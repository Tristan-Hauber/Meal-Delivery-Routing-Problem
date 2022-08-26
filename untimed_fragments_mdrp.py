#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 09:24:49 2022

@author: Tristan Hauber
"""

from __future__ import annotations

from gurobi import Model, quicksum, GRB
from classes import Group, Arc, Data, Order, UntimedFragmentPath, Courier, Fragment
from typing import List


class UntimedFragmentsMDRP(Model):
    """
    A meal delivery routing problem model, solved using untimed path fragments.

    Attributes
    ----------
    group : Group
        The courier group on which the subproblem is defined.
    couriers : List[Courier]
        The couriers in the group.
    arcs : List[Arc]
        The arcs used to solve the model.
    orders : List[Order]
        The orders needed to be delivered.

    """

    _model_from_group_and_orders = dict()

    def __init__(self, group: Group, arcs: List[Arc], orders: List[Order], cost_penalty: int = 10000,
                 cost_penalty_active: bool = True, time_limit: int = 5, gap: int = 9989, save_model: bool = False,
                 *args: object, **kwargs: object) -> None:
        """
        Create a new meal delivery routing problem model.

        Parameters
        ----------
        group : Group
            The courier group to deliver orders.
        arcs : List[Arc]
            The list of possible arcs used to deliver orders.
        orders : List[Order]
            The list of orders needing to be delivered.

        Returns
        -------
        None.

        """
        super().__init__(*args, **kwargs)
        """ ========== SETUP ========== """
        if save_model:
            UntimedFragmentsMDRP._model_from_group_and_orders[(group, frozenset(orders))] = self
        # Model
        self.setParam('OutputFlag', 0)
        self.setParam('TimeLimit', time_limit)
        self.setParam('MIPGapAbs', gap)
        """ ========== SETS ========== """
        self._group = group
        self._couriers = group.couriers
        self._arcs = arcs
        self._orders = orders
        """ ========== VARIABLES ========== """
        # Payment to each courier
        self._courier_payments = {courier: self.addVar() for courier in self._couriers}
        # Time of departure for each arc
        self._timings = {arc: self.addVar() for arc in self._arcs}
        # If the order is delivered
        self._deliveries = {
            order: self.addVar(vtype=GRB.BINARY) for order in self._orders
        }
        # If the arc is serviced
        self._serviced = {arc: self.addVar(vtype=GRB.BINARY) for arc in self._arcs}
        # Courier assigned to arc
        self._assignments = {
            (courier, arc): self.addVar(vtype=GRB.BINARY)
            for courier in self._couriers
            for arc in self._arcs
        }
        # Arc2 follows arc1
        self._successors = {
            (arc1, arc2): self.addVar(vtype=GRB.BINARY)
            for arc1 in arcs
            for arc2 in arcs
            if arc2 != arc1
            and arc1.arrival_location == arc2.departure_location
            and arc1.earliest_departure_time + arc1.travel_time <= arc2.latest_departure_time
        }

        """ ========== OBJECTIVE ========== """
        self.setObjective(
            quicksum(self._courier_payments[courier] for courier in self._couriers)
            + quicksum(
                (1 - self._deliveries[order]) * cost_penalty for order in self._orders
            )
        )

        """ ========== CONSTRAINTS ========== """
        # Pay couriers for their time
        self._time_payments = {
            courier: self.addConstr(
                self._courier_payments[courier]
                >= (courier.off_time - courier.on_time) / 60 * Data.MIN_PAY_PER_HOUR
            )
            for courier in self._couriers
        }
        # Pay couriers for their deliveries
        self._delivery_payments = {
            courier: self.addConstr(
                self._courier_payments[courier]
                >= quicksum(
                    self._assignments[courier, arc]
                    * len(arc.order_list)
                    * Data.PAY_PER_DELIVERY
                    for arc in self._arcs
                )
            )
            for courier in self._couriers
        }
        # All orders must be delivered
        if cost_penalty_active is False:
            self._deliver_orders = {
                order: self.addConstr(self._deliveries[order] == 1)
                for order in self._orders
            }
        # If delivered, an arc must service the order
        self._delivered = {
            order: self.addConstr(
                self._deliveries[order]
                == quicksum(
                    self._serviced[arc] for arc in self._arcs if order in arc.order_list
                )
            )
            for order in self._orders
        }
        # All activated arcs are assigned to a courier
        self._arcs_covered = {
            arc: self.addConstr(
                self._serviced[arc]
                == quicksum(
                    self._assignments[courier, arc]
                    for courier in self._couriers
                    if (courier, arc) in self._assignments
                )
            )
            for arc in self._arcs
        }
        # Predecessor and successor setup
        predecessors = set()
        successors = set()
        for (arc1, arc2) in self._successors:
            predecessors.add(arc1)
            successors.add(arc2)
        # All activated non-exit arcs must have a successor
        self._has_successor = {
            arc1: self.addConstr(
                self._serviced[arc1]
                == quicksum(
                    self._successors[arc1, arc2]
                    for arc2 in arcs
                    if (arc1, arc2) in self._successors
                )
            )
            for arc1 in predecessors
        }
        # All activated non-entry arcs must have a predecessor
        self._has_predecessor = {
            arc2: self.addConstr(
                self._serviced[arc2] == quicksum(
                    self._successors[arc1, arc2] for arc1 in arcs if (arc1, arc2) in self._successors)) for arc2 in
            successors}
        # All successors must start late enough for the previous to finish
        self._successor_timings = {
            (arc1, arc2): self.addConstr(
                self._timings[arc1] + arc1.travel_time
                <= self._timings[arc2]
                + (1 - self._successors[arc1, arc2])
                * (
                        arc1.latest_departure_time
                        + arc1.travel_time
                        - arc2.earliest_departure_time
                )
            )
            for (arc1, arc2) in self._successors
        }
        # If a courier services an arc, it must also service its successor
        self._courier_service_successor = {
            (arc1, arc2, courier): self.addConstr(
                self._assignments[courier, arc1] + self._successors[arc1, arc2] - 1
                <= self._assignments[courier, arc2]
            )
            for courier in self._couriers
            for (arc1, arc2) in self._successors
        }
        # If an arc is serviced, it must have a courier servicing it
        self._servicing_requires_courier = {
            arc: self.addConstr(
                self._serviced[arc]
                == quicksum(self._assignments[courier, arc] for courier in self._couriers)
            )
            for arc in self._arcs
        }
        # An arc must be serviced after it is ready
        self._begin_on_time = {
            arc: self.addConstr(self._timings[arc] >= arc.earliest_departure_time)
            for arc in self._arcs
        }
        # An arc must be serviced before it is too late
        self._begin_in_time = {
            arc: self.addConstr(self._timings[arc] <= arc.latest_departure_time)
            for arc in self._arcs
        }
        # A courier may not service more than one entry arc
        self._start_once = {courier: self.addConstr(
            quicksum(
                self._serviced[arc]
                for arc in self._serviced
                if arc.departure_location == courier)
            <= 1)
            for courier in self._couriers}
        """ ========== Attributes ========== """
        self._paths = None

    def get_arc_network_paths(self) -> List[UntimedFragmentPath]:
        """
        Create a list of paths through the network.

        Each path through the network is of type UntimedFragmentPath, and this
        list of paths partitions the entire arc space.

        This function finds all activated untimed path fragments and transitions
        between those fragments. It chooses an arbitrary starting fragment, then
        follows the transitions until it comes to an exit fragment.

        At the end of the function, it should have covered all activated
        fragments and transitions.

        Returns
        -------
        List[UntimedFragmentPath]
            A spanning list of UntimedFragmentPath through the network.

        """
        if self._paths is not None:
            return self._paths
        activated_untimed_path_fragments = set()
        activated_starting_untimed_path_fragments = set()
        for untimed_path_fragment in self._arcs:
            if self._serviced[untimed_path_fragment].x > 0.9:
                if type(untimed_path_fragment.departure_location) == Courier:
                    activated_starting_untimed_path_fragments.add(untimed_path_fragment)
                else:
                    activated_untimed_path_fragments.add(untimed_path_fragment)
        activated_transitions = set()
        for transition in self._successors:
            if self._successors[transition].x > 0.9:
                activated_transitions.add(transition)
        self._paths = list()
        # Get a path start
        while len(activated_starting_untimed_path_fragments) > 0:
            path = list()
            path.append(activated_starting_untimed_path_fragments.pop())
            # Follow it to the end
            while True:
                found_following_fragment = False
                for transition in activated_transitions:
                    if transition[0] == path[-1]:
                        new_untimed_path_fragment = transition[1]
                        path.append(new_untimed_path_fragment)
                        activated_untimed_path_fragments.remove(new_untimed_path_fragment)
                        activated_transitions.remove(transition)
                        found_following_fragment = True
                        break
                if type(path[-1].arrival_location) == Group:
                    break
                assert found_following_fragment
            untimed_fragment_path = UntimedFragmentPath(path)
            self._paths.append(untimed_fragment_path)
        assert len(activated_transitions) == 0
        assert len(activated_untimed_path_fragments) == 0
        return self._paths

    def convert_to_timed_path_fragments(self) -> List[Fragment]:
        """Get a list of timed path fragments from the network."""
        if self._paths is None:
            self.get_arc_network_paths()
        return Fragment.get_timed_fragments_from_untimed_fragment_network(self._paths)

    def get_undelivered_orders(self) -> List[Order]:
        """Return a list of all undelivered orders in the model."""
        self.update()
        return list(order for order in self._deliveries if self._deliveries[order].x < 0.1)

    @staticmethod
    def get_arc_model(group: Group, orders: List[Order]) -> UntimedFragmentsMDRP:
        """Return the model built on the given courier group and order set."""
        if (group, frozenset(orders)) in UntimedFragmentsMDRP._model_from_group_and_orders:
            return UntimedFragmentsMDRP._model_from_group_and_orders[(group, frozenset(orders))]

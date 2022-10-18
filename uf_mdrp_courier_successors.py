#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 09:24:49 2022

@author: Tristan Hauber
"""

from __future__ import annotations

from gurobipy import Model, quicksum, GRB
from classes import Group, Arc, Data, Order, UntimedFragmentPath, Courier, Fragment
from typing import List, Set, Dict, Tuple

import time
import math
import itertools


class UFMDRP2(Model):
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
                 cost_penalty_active: bool = False, time_limit: int = GRB.INFINITY, gap: int = 1e-10,
                 save_model: bool = False, output_type: str = "Summary", *args: object, **kwargs: object) -> None:
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
        self._model_initiation = time.time()
        super().__init__(*args, **kwargs)
        """ ========== SETUP ========== """
        if save_model:
            UFMDRP2._model_from_group_and_orders[(group, frozenset(orders))] = self
        # Model
        self.setParam('OutputFlag', 0)
        self.setParam('TimeLimit', time_limit)
        self.setParam('MIPGapAbs', gap)
        """ ========== SETS ========== """
        self._group = group
        self._couriers = group.couriers
        self._arcs = arcs
        self._orders = orders
        """ ========== OPTIMISATION ========== """
        self._couriers_for_arc: Dict[Arc: Set[Courier]] = {arc: set() for arc in self._arcs}
        self._arcs_for_courier: Dict[Courier: Set[Arc]] = {courier: set() for courier in self._couriers}
        self._entry_arcs_for_courier: Dict[Courier: Set[Arc]] = {courier: set() for courier in self._couriers}
        self._arcs_for_order: Dict[Order: Set[Arc]] = {order: set() for order in self._orders}
        self._successors_of_arc: Dict[Arc: Set[Arc]] = {arc: set() for arc in self._arcs if
                                                        type(arc.arrival_location) != Group}
        self._predecessors_of_arc: Dict[Arc: Set[Arc]] = {arc: set() for arc in self._arcs if
                                                          type(arc.departure_location) != Courier}
        self._entry_arcs: Set[Arc] = set()
        self._exit_arcs: Set[Arc] = set()
        # Find possible servicing couriers for the arc
        for arc in self._arcs:
            if type(arc.departure_location) == Courier:
                self._couriers_for_arc[arc].add(arc.departure_location)
                self._arcs_for_courier[arc.departure_location].add(arc)
                self._entry_arcs_for_courier[arc.departure_location].add(arc)
                self._entry_arcs.add(arc)
            else:
                for courier in self._couriers:
                    if courier.get_earliest_arrival_at(arc.departure_location) <= arc.latest_departure_time:
                        self._couriers_for_arc[arc].add(courier)
                        self._arcs_for_courier[courier].add(arc)

                if type(arc.arrival_location) == Group:
                    self._exit_arcs.add(arc)

            # Find orders serviced by the arc
            for order in arc.order_list:
                self._arcs_for_order[order].add(arc)

        # Find possible pred-succ pairs
        for (arc1, arc2) in itertools.combinations(self._arcs, 2):
            if set(arc1.order_list).intersection(arc2.order_list) == set():
                if arc1.arrival_location == arc2.departure_location \
                        and arc1.earliest_departure_time + arc1.travel_time <= arc2.latest_departure_time:
                    self._successors_of_arc[arc1].add(arc2)
                    self._predecessors_of_arc[arc2].add(arc1)
                if arc2.arrival_location == arc1.departure_location \
                        and arc2.earliest_departure_time + arc2.travel_time <= arc1.latest_departure_time:
                    self._successors_of_arc[arc2].add(arc1)
                    self._predecessors_of_arc[arc1].add(arc2)
        """ ========== VARIABLES ========== """
        # Payment to each courier
        self._courier_payments = {courier: self.addVar() for courier in self._couriers}
        # Time of departure for each arc
        self._timings = {arc: self.addVar() for arc in self._arcs}
        # If the order is delivered
        self._deliveries = {
            order: self.addVar(vtype=GRB.BINARY) for order in self._orders
        }
        # Courier assigned to arc
        self._assignments = {
            (courier, arc): self.addVar(vtype=GRB.BINARY)
            for courier in self._couriers
            for arc in self._arcs_for_courier[courier]
        }
        # Arc2 follows arc1
        self._successors = {
            (courier, arc1, arc2): self.addVar(vtype=GRB.BINARY)
            for arc1 in self._successors_of_arc
            for arc2 in self._successors_of_arc[arc1]
            for courier in self._couriers_for_arc[arc1]
            if courier in self._couriers_for_arc[arc2]
        }
        if output_type == "Model Size":
            self.update()
            print(
                f'{self.getAttr(GRB.Attr.NumVars)} variables created, t={math.ceil(time.time() - self._model_initiation)}')

        """ ========== OBJECTIVE ========== """
        self.setObjective(
            quicksum(self._courier_payments[courier] for courier in self._couriers)
            + quicksum(
                (1 - self._deliveries[order]) * cost_penalty for order in self._orders
            )
        )

        """ ========== CONSTRAINTS ========== """
        # Pay couriers for their deliveries
        self._delivery_payments = {
            courier: self.addConstr(
                self._courier_payments[courier]
                >= quicksum(
                    self._assignments[courier, arc]
                    * len(arc.order_list)
                    * Data.PAY_PER_DELIVERY
                    for arc in self._arcs_for_courier[courier]
                )
            )
            for courier in self._couriers
        }
        # Pay couriers for their time
        self._time_payments = {
            courier: self.addConstr(
                self._courier_payments[courier]
                >= (courier.off_time - courier.on_time) / 60 * Data.MIN_PAY_PER_HOUR
            )
            for courier in self._couriers
        }
        # All orders must be delivered
        if cost_penalty_active is False:
            self._deliver_orders = {
                order: self.addConstr(self._deliveries[order] == 1)
                for order in self._orders
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
        # If delivered, an arc must service the order
        self._activate_arc_for_order = {
            order: self.addConstr(quicksum(self._assignments[courier, arc]
                                           for (courier, arc) in self._assignments
                                           if order in arc.order_list)
                                  == self._deliveries[order])
            for order in self._orders}
        # Each activated arc must have a successor by the same courier
        self._have_successor_constraint = {
            (courier, arc1): self.addConstr(
                self._assignments[(courier, arc1)]
                == quicksum(self._successors[(courier, arc1, arc2)]
                            for arc2 in self._successors_of_arc[arc1]
                            if (courier, arc1, arc2) in self._successors))
            for (courier, arc1) in self._assignments
            if type(arc1.arrival_location) != Group}
        # Each activated arc must have a predecessor by the same courier
        self._have_predecessor_constraint = {
            (courier, arc2): self.addConstr(
                self._assignments[courier, arc2]
                == quicksum(self._successors[courier, arc1, arc2]
                            for arc1 in self._predecessors_of_arc[arc2]
                            if (courier, arc1, arc2) in self._successors))
            for (courier, arc2) in self._assignments
            if type(arc2.departure_location) != Courier}
        # All successors must start late enough for the previous to finish
        self._successor_timings = {
            (arc1, arc2): self.addConstr(
                self._timings[arc1] + arc1.travel_time
                <= self._timings[arc2]
                + (1 - quicksum(self._successors[courier, arc1, arc2]
                                for courier in self._couriers
                                if (courier, arc1, arc2) in self._successors))
                * (
                        arc1.latest_departure_time
                        + arc1.travel_time
                        - arc2.earliest_departure_time
                )
            )
            for arc1 in self._successors_of_arc
            for arc2 in self._successors_of_arc[arc1]
        }
        # A courier may not service more than one entry arc
        self._start_once = {courier: self.addConstr(
            quicksum(
                self._assignments[courier, arc]
                for arc in self._entry_arcs_for_courier[courier])
            <= 1)
            for courier in self._couriers}
        if output_type == "Model Size":
            self.update()
            print(
                f'{self.getAttr(GRB.Attr.NumConstrs)} constraints created, t={math.ceil(time.time() - self._model_initiation)}')
        """ ========== Attributes ========== """
        self._paths = None

    def get_arc_network_paths(self) -> List[UntimedFragmentPath]:
        """
        Create a list of paths through the network.

        Each path through the network is of type UntimedFragmentPath, and this
        list of paths partitions the entire arc space.

        This function works by finding all activated transitions in the
        self._successors variable, and sorts them by courier. Then for each
        courier, the transitions are followed, starting at an entry arc, and
        finishing at an exit arc.

        Checks are run to ensure all arcs are covered.

        Returns
        -------
        List[UntimedFragmentPath]
            A spanning list of UntimedFragmentPath through the network.

        """
        # If already calculated, no need to do it again
        if self._paths is not None:
            return self._paths

        # Find all activated transitions, and sort them by courier
        activated_transitions_by_courier: Dict[Courier: Set[Tuple[Arc, Arc]]] = {courier: set() for courier in self._couriers}
        for (courier, arc1, arc2) in self._successors:
            if self._successors[courier, arc1, arc2].x > 0.9:
                activated_transitions_by_courier[courier].add((arc1, arc2))

        # Courier by courier, follow the transitions from entry to exit, add path
        self._paths: List[UntimedFragmentPath] = list()
        for courier in self._couriers:
            path: List[Arc] = list()
            working_transitions = activated_transitions_by_courier[courier]
            # Initiate start arc as the entry arc
            for arc in self._arcs_for_courier[courier]:
                if self._assignments[courier, arc].x > 0.9 and arc.departure_location == courier:
                    path.append(arc)
            # Check that exactly one entry arc was found
            assert len(path) == 1 or len(working_transitions) == 0
            # Find the next transition
            while len(working_transitions) > 0:
                for (arc1, arc2) in working_transitions:
                    # Check that the transition starts with the last arc
                    if arc1 == path[-1]:
                        path.append(arc2)
                        working_transitions.remove((arc1, arc2))
                        break
                else:
                    # No valid transitions found
                    assert False
            # Only add the path if it is non-empty and valid
            if len(path) > 0:
                # The path should end at the group depot
                assert type(path[-1].arrival_location) == Group
                # Add the path to the list of paths
                self._paths.append(UntimedFragmentPath(path))

        # Return completed list of paths
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
    def get_arc_model(group: Group, orders: List[Order]) -> UFMDRP2:
        """Return the model built on the given courier group and order set."""
        if (group, frozenset(orders)) in UFMDRP2._model_from_group_and_orders:
            return UFMDRP2._model_from_group_and_orders[(group, frozenset(orders))]

    def print_infeasible_constraints(self) -> None:
        """Print the name of all IIS constraints."""
        self.computeIIS()
        print("\nIIS Constraints include:")
        for courier in self._delivery_payments:
            if self._delivery_payments[courier].IISConstr:
                print("Delivery payment constraints")
                break
        for courier in self._time_payments:
            if self._time_payments[courier].IISConstr:
                print("Time payment constraints")
                break
        for order in self._deliver_orders:
            if self._deliver_orders[order].IISConstr:
                print("Order delivery constraints")
                break
        for (arc1, arc2) in self._successor_timings:
            if self._successor_timings[(arc1, arc2)].IISConstr:
                print("Arc successor timing constraints")
                break
        for arc in self._begin_on_time:
            if self._begin_on_time[arc].IISConstr:
                print("Arc earliest departure time constraints")
                break
        for arc in self._begin_in_time:
            if self._begin_in_time[arc].IISConstr:
                print("Arc latest departure time constraints")
                break
        for courier in self._start_once:
            if self._start_once[courier].IISConstr:
                print("Couriers start once constraints")
                break
        for order in self._activate_arc_for_order:
            if self._activate_arc_for_order[order].IISConstr:
                print("Activation of arcs for order delivery constraints")
                break
        for (courier, arc1) in self._have_successor_constraint:
            if self._have_successor_constraint[(courier, arc1)].IISConstr:
                print("Requirement for an arc to have a successor")
                break
        for (courier, arc2) in self._have_predecessor_constraint:
            if self._have_predecessor_constraint[(courier, arc2)].IISConstr:
                print("Requirement for an arc to have a predecessor")
                break
        print()

    def optimize(self):
        super().optimize()
        # print(f'Optimisation complete, t={math.ceil(time.time() - self._model_initiation)}')
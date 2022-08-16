#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 09:24:49 2022

@author: Tristan Hauber
"""

from gurobi import Model, quicksum, GRB
from classes import Group, Arc, Data, Order, UntimedFragmentPath, Courier, Fragment
from typing import List


class UntimedFragmentsMDRP():
    """
    A meal delivery routing problem model, solved using untimed path fragments.

    Attributes
    ----------
    couriers : List[Courier]
        The couriers in the group.
    arcs : List[Arc]
        The arcs used to solve the model.
    orders : List[Order]
        The orders needed to be delivered.

    """

    def __init__(
        self,
        group: Group,
        arcs: List[Arc],
        orders: List[Order],
        cost_penalty=10000,
        cost_penalty_active=True,
    ):
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
        self.uf_mdrp = Model('Untimed Fragments MDRP')
        self.uf_mdrp.setParam('OutputFlag', 0)
        """ ========== SETS ========== """
        self.couriers = group.couriers
        self.arcs = arcs
        self.orders = orders

        """ ========== VARIABLES ========== """
        # Payment to each courier
        self.courier_payments = {courier: self.uf_mdrp.addVar() for courier in self.couriers}
        # Time of departure for each arc
        self.timings = {arc: self.uf_mdrp.addVar() for arc in self.arcs}
        # If the order is delivered
        self.deliveries = {
            order: self.uf_mdrp.addVar(vtype=GRB.BINARY) for order in self.orders
        }
        # If the arc is serviced
        self.serviced = {arc: self.uf_mdrp.addVar(vtype=GRB.BINARY) for arc in self.arcs}
        # Courier assigned to arc
        self.assignments = {
            (courier, arc): self.uf_mdrp.addVar(vtype=GRB.BINARY)
            for courier in self.couriers
            for arc in self.arcs
        }
        # Arc2 follows arc1
        self.successors = {
            (arc1, arc2): self.uf_mdrp.addVar(vtype=GRB.BINARY)
            for arc1 in arcs
            for arc2 in arcs
            if arc2 != arc1
            and arc1.arrival_location == arc2.departure_location
            and arc1.earliest_departure_time + arc1.travel_time
            <= arc2.latest_departure_time
        }

        """ ========== OBJECTIVE ========== """
        self.uf_mdrp.setObjective(
            quicksum(self.courier_payments[courier] for courier in self.couriers)
            + quicksum(
                (1 - self.deliveries[order]) * cost_penalty for order in self.orders
            )
        )

        """ ========== CONSTRAINTS ========== """
        # Pay couriers for their time
        self.time_payments = {
            courier: self.uf_mdrp.addConstr(
                self.courier_payments[courier]
                >= (courier.off_time - courier.on_time) / 60 * Data.MIN_PAY_PER_HOUR
            )
            for courier in self.couriers
        }
        # Pay couriers for their deliveries
        self.delivery_payments = {
            courier: self.uf_mdrp.addConstr(
                self.courier_payments[courier]
                >= quicksum(
                    self.assignments[courier, arc]
                    * len(arc.order_list)
                    * Data.PAY_PER_DELIVERY
                    for arc in self.arcs
                )
            )
            for courier in self.couriers
        }
        # All orders must be delivered
        if cost_penalty_active is False:
            self.deliver_orders = {
                order: self.uf_mdrp.addConstr(self.deliveries[order] == 1)
                for order in self.orders
            }
        # If delivered, an arc must service the order
        self.delivered = {
            order: self.uf_mdrp.addConstr(
                self.deliveries[order]
                == quicksum(
                    self.serviced[arc] for arc in self.arcs if order in arc.order_list
                )
            )
            for order in self.orders
        }
        # All activated arcs are assigned to a courier
        self.arcs_covered = {
            arc: self.uf_mdrp.addConstr(
                self.serviced[arc]
                == quicksum(
                    self.assignments[courier, arc]
                    for courier in self.couriers
                    if (courier, arc) in self.assignments
                )
            )
            for arc in self.arcs
        }
        # Predecessor and successor setup
        predecessors = set()
        successors = set()
        for (arc1, arc2) in self.successors:
            predecessors.add(arc1)
            successors.add(arc2)
        # All activated non-exit arcs must have a successor
        self.has_successor = {
            arc1: self.uf_mdrp.addConstr(
                self.serviced[arc1]
                == quicksum(
                    self.successors[arc1, arc2]
                    for arc2 in arcs
                    if (arc1, arc2) in self.successors
                )
            )
            for arc1 in predecessors
        }
        # All activated non-entry arcs must have a predecessor
        self.has_predecessor = {arc2: self.uf_mdrp.addConstr(self.serviced[arc2] == quicksum(self.successors[arc1, arc2] for arc1 in arcs if (arc1, arc2) in self.successors)) for arc2 in successors}
        # All successors must start late enough for the previous to finish
        self.succ_timings = {
            (arc1, arc2): self.uf_mdrp.addConstr(
                self.timings[arc1] + arc1.travel_time
                <= self.timings[arc2]
                + (1 - self.successors[arc1, arc2])
                * (
                    arc1.latest_departure_time
                    + arc1.travel_time
                    - arc2.earliest_departure_time
                )
            )
            for (arc1, arc2) in self.successors
        }
        # If a courier services an arc, it must also service its successor
        self.cour_serv_succ = {
            (arc1, arc2, courier): self.uf_mdrp.addConstr(
                self.assignments[courier, arc1] + self.successors[arc1, arc2] - 1
                <= self.assignments[courier, arc2]
            )
            for courier in self.couriers
            for (arc1, arc2) in self.successors
        }
        # If an arc is serviced, it must have a courier servicing it
        self.serv_requ_cour = {
            arc: self.uf_mdrp.addConstr(
                self.serviced[arc]
                == quicksum(self.assignments[courier, arc] for courier in self.couriers)
            )
            for arc in self.arcs
        }
        # An arc must be serviced after it is ready
        self.begin_on_time = {
            arc: self.uf_mdrp.addConstr(self.timings[arc] >= arc.earliest_departure_time)
            for arc in self.arcs
        }
        # An arc must be serviced before it is too late
        self.begin_in_time = {
            arc: self.uf_mdrp.addConstr(self.timings[arc] <= arc.latest_departure_time)
            for arc in self.arcs
        }
        """ ========== Attributes ========== """
        self.paths = None

    def get_untimed_path_fragment_path(self) -> List[UntimedFragmentPath]:
        """
        Create a list of paths through the network.

        Each path through the network is of type UntimedFragmentPath, and this list of paths is spanning.

        This function finds all activated untimed path fragments and transitions between those fragments. It chooses an arbitrary starting fragment, then follows the transitions until it comes to an exit fragment.

        At the end of the function, it should have covered all activated fragments and transitions.

        Returns
        -------
        List[UntimedFragmentPath]
            A spanning list of UntimedFragmentPath through the network.

        """
        if self.paths is not None:
            return self.paths
        activated_untimed_path_fragments = set()
        activated_starting_untimed_path_fragments = set()
        for untimed_path_fragment in self.arcs:
            if self.serviced[untimed_path_fragment].x > 0.9:
                if type(untimed_path_fragment.departure_location) == Courier:
                    activated_starting_untimed_path_fragments.add(untimed_path_fragment)
                else:
                    activated_untimed_path_fragments.add(untimed_path_fragment)
        activated_transitions = set()
        for transition in self.successors:
            if self.successors[transition].x > 0.9:
                activated_transitions.add(transition)
        self.paths = list()
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
            self.paths.append(untimed_fragment_path)
        assert len(activated_transitions) == 0
        assert len(activated_untimed_path_fragments) == 0
        return self.paths

    def convert_to_timed_path_fragments(self) -> List[Fragment]:
        """Get a list of timed path fragments from the network."""
        if self.paths is None:
            self.get_untimed_path_fragment_path()
        return Fragment.get_timed_fragments_from_untimed_fragment_network(self.paths)

    def optimize(self, callback=None):
        """Optimize the gurobi model, with optional callback possibilities."""
        self.uf_mdrp.optimize(callback)
        self.objVal = self.uf_mdrp.objVal

    def get_undelivered_orders(self) -> List[Order]:
        """Return a list of all undelivered orders in the model."""
        self.uf_mdrp.update()
        return list(order for order in self.deliveries if self.deliveries[order].x < 0.1)

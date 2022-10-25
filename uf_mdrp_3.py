#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 09:24:49 2022

@author: Tristan Hauber
"""

from __future__ import annotations

from classes import Group, Arc, Order, Courier, Data, UntimedFragmentPath, Fragment
from typing import List, Dict, Set, Tuple
from gurobipy import Model, GRB, quicksum

import itertools
import math
import time


class UFMDRP3(Model):
    """A courier group sub-model, solved using a callback."""

    _model_from_group_and_orders = dict()

    def __init__(self, group: Group, arcs: List[Arc], orders: List[Order], cost_penalty: int = 10000,
                 cost_penalty_active: bool = False, time_limit: int = GRB.INFINITY, gap: int = 1e-10,
                 save_model: bool = False, output_type: str = "None", *args: object, **kwargs: object) -> None:
        """Create a new MDRP group submodel."""
        super().__init__(*args, **kwargs)
        """ ========== SETUP ========== """
        if save_model:
            UFMDRP3._model_from_group_and_orders[(group, frozenset(orders))] = self
        # Model
        if output_type != "Full":
            self.setParam('OutputFlag', 0)
        self.setParam('TimeLimit', time_limit)
        self.setParam('MIPGapAbs', gap)
        self.setParam('LazyConstraints', 1)
        self._model_initiation = time.time()
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
        # print(self._couriers)
        # print(self._arcs)
        # print(self._orders)
        self._courier_payment_variables = {courier:
                                           self.addVar(lb=0, ub=GRB.INFINITY, obj=0, vtype=GRB.CONTINUOUS,
                                                       name=f'{courier} payment variable', column=None)
                                           for courier in self._couriers}
        self._order_delivery_variables = {order:
                                          self.addVar(lb=0, ub=1, obj=0, vtype=GRB.BINARY,
                                                      name=f'{order} delivery variable', column=None)
                                          for order in self._orders}
        self._courier_assignment_variables = {(courier, arc):
                                              self.addVar(lb=0, ub=1, obj=0, vtype=GRB.BINARY,
                                                          name=f'{(courier, arc)} assignment variable', column=None)
                                              for courier in self._couriers
                                              for arc in self._arcs_for_courier[courier]}
        self._courier_successor_variables = {(courier, arc1, arc2):
                                             self.addVar(lb=0, ub=1, obj=0, vtype=GRB.BINARY,
                                                         name=f'{(courier, arc1, arc2)} successor variable',
                                                         column=None)
                                             for arc1 in self._successors_of_arc
                                             for arc2 in self._successors_of_arc[arc1]
                                             for courier in self._couriers_for_arc[arc1]
                                             if courier in self._couriers_for_arc[arc2]}
        if output_type != "None":
            self.update()
            print(
                f'{self.getAttr(GRB.Attr.NumVars)} variables created, '
                f't={math.ceil(time.time() - self._model_initiation)}')

        """ ========== OBJECTIVE ========== """
        self.setObjective(
            quicksum(self._courier_payment_variables[courier] for courier in self._couriers)
            + quicksum(
                (1 - self._order_delivery_variables[order]) * cost_penalty for order in self._orders
            )
        )

        """ ========== CONSTRAINTS ========== """
        # Pay couriers for their deliveries
        self._delivery_payments = {
            courier: self.addConstr(
                self._courier_payment_variables[courier]
                >= quicksum(
                    self._courier_assignment_variables[courier, arc]
                    * len(arc.order_list)
                    * Data.PAY_PER_DELIVERY
                    for arc in self._arcs_for_courier[courier]
                ), name=f'{courier} delivery payment constraint'
            )
            for courier in self._couriers
        }
        # Pay couriers for their time
        self._time_payments = {
            courier: self.addConstr(
                self._courier_payment_variables[courier]
                >= (courier.off_time - courier.on_time) / 60 * Data.MIN_PAY_PER_HOUR,
                name=f'{courier} time payment constraint')
            for courier in self._couriers
        }
        # All orders must be delivered
        if cost_penalty_active is False:
            self._deliver_orders = {
                order: self.addConstr(self._order_delivery_variables[order] == 1, f'{order} delivery constraint')
                for order in self._orders
            }
        # If delivered, an arc must service the order
        self._activate_arc_for_order = {
            order: self.addConstr(quicksum(self._courier_assignment_variables[courier, arc]
                                           for (courier, arc) in self._courier_assignment_variables
                                           if order in arc.order_list)
                                  == self._order_delivery_variables[order], name=f'{order} assignment constraint')
            for order in self._orders}
        # Each activated arc must have a successor by the same courier
        self._have_successor_constraint = {
            (courier, arc1): self.addConstr(
                self._courier_assignment_variables[(courier, arc1)]
                == quicksum(self._courier_successor_variables[(courier, arc1, arc2)]
                            for arc2 in self._successors_of_arc[arc1]
                            if (courier, arc1, arc2) in self._courier_successor_variables),
                name=f'{courier, arc1} successor constraint')
            for (courier, arc1) in self._courier_assignment_variables
            if type(arc1.arrival_location) != Group}
        # Each activated arc must have a predecessor by the same courier
        self._have_predecessor_constraint = {
            (courier, arc2): self.addConstr(
                self._courier_assignment_variables[courier, arc2]
                == quicksum(self._courier_successor_variables[courier, arc1, arc2]
                            for arc1 in self._predecessors_of_arc[arc2]
                            if (courier, arc1, arc2) in self._courier_successor_variables),
                name=f'{courier, arc2} predecessor constraint')
            for (courier, arc2) in self._courier_assignment_variables
            if type(arc2.departure_location) != Courier}
        # A courier may not service more than one entry arc
        self._start_once = {courier: self.addConstr(
            quicksum(
                self._courier_assignment_variables[courier, arc]
                for arc in self._entry_arcs_for_courier[courier])
            <= 1, name=f'{courier} start count constraint')
            for courier in self._couriers}
        if output_type != "None":
            self.update()
            print(
                f'{self.getAttr(GRB.Attr.NumConstrs)} constraints created, '
                f't={math.ceil(time.time() - self._model_initiation)}')
        """ ========== Attributes ========== """
        self._paths = None

    def get_arc_network_paths(self) -> List[UntimedFragmentPath]:
        """
        Create a list of paths through the network.

        Each path through the network is of type UntimedFragmentPath, and this
        list of paths partitions the entire arc space.

        This function works by finding all activated transitions in the
        self._courier_successor_variables variable, and sorts them by courier. Then for each
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
        activated_transitions_by_courier: Dict[Courier: Set[Tuple[Arc, Arc]]] = {courier: set() for courier in
                                                                                 self._couriers}
        for (courier, arc1, arc2) in self._courier_successor_variables:
            if self._courier_successor_variables[courier, arc1, arc2].x > 0.9:
                activated_transitions_by_courier[courier].add((arc1, arc2))

        # Courier by courier, follow the transitions from entry to exit, add path
        self._paths: List[UntimedFragmentPath] = list()
        for courier in self._couriers:
            path: List[Arc] = list()
            working_transitions = activated_transitions_by_courier[courier]
            # Initiate start arc as the entry arc
            for arc in self._arcs_for_courier[courier]:
                if self._courier_assignment_variables[courier, arc].x > 0.9 and arc.departure_location == courier:
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
        return list(order for order in self._order_delivery_variables if self._order_delivery_variables[order].x < 0.1)

    @staticmethod
    def get_arc_model(group: Group, orders: List[Order]) -> UFMDRP3:
        """Return the model built on the given courier group and order set."""
        if (group, frozenset(orders)) in UFMDRP3._model_from_group_and_orders:
            return UFMDRP3._model_from_group_and_orders[(group, frozenset(orders))]

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

    def optimize(self) -> None:
        """Optimise the submodel, using a callback."""

        def trace_back_arcs(transitions: Set[Tuple[Arc, Arc]], starting_arc: Arc) -> Set[Tuple[Arc, Arc]]:
            """
            Follow the transitions backward to find a minimal infeasible path.

            This function returns two values in a tuple, the first is the set
            of transitions defining the minimal infeasible path. The second is
            the second arc in the path, from the start.
            """
            current_arc = starting_arc
            current_time = current_arc.latest_departure_time
            infeasible_transitions: Set[Tuple[Arc, Arc]] = set()
            # Trace backwards, updating time
            while True:
                print("Hi")
                for (arc1, arc2) in transitions:
                    if arc2 == current_arc:
                        infeasible_transitions.add((arc1, arc2))
                        current_time -= arc1.travel_time
                        if current_time < arc1.earliest_departure_time:
                            # Minimal infeasible path found
                            return infeasible_transitions
                        else:
                            # Valid transition, move backwards along path
                            current_arc = arc1

                        # Transition found, move to next transition
                        break

                else:
                    # No valid transition found
                    assert False

        def callback(model: Model, where: int) -> None:
            """Verify that the found given path is valid."""
            """
            We solely have to assign a time to every arc. We do this by starting
            with the first arc in the path, then follow the succession variables,
            assigning times as we go. If we find that the timing is impossible,
            remove the subset of arcs as infeasible.
            """
            if where == GRB.Callback.MIPSOL:
                # Find all activated transitions
                activated_transitions: Dict[Courier, Set[Tuple[Arc, Arc]]] = {courier: set() for courier in
                                                                              self._couriers}
                for (courier, arc1, arc2) in self._courier_successor_variables:
                    if self.cbGetSolution(self._courier_successor_variables[courier, arc1, arc2]) > 0.9:
                        activated_transitions[courier].add((arc1, arc2))

                # Iterate through each courier path
                for courier in self._couriers:

                    # Retrieve activated transitions for the courier
                    courier_transitions = activated_transitions[courier]
                    remaining_transitions = set(activated_transitions[courier])

                    if len(courier_transitions) == 0:
                        continue

                    # Find the entry arc
                    current_arc: Arc
                    current_time: int
                    for arc in self._entry_arcs_for_courier[courier]:
                        if self.cbGetSolution(self._courier_assignment_variables[courier, arc]) > 0.9:
                            current_arc = arc
                            current_time = current_arc.earliest_departure_time
                            break
                    else:
                        # No initial arc found, but arcs exist
                        # Cut off set of transitions
                        # print('No entry arc')
                        for courier2 in self._couriers:
                            if all((courier2, arc1, arc2) in self._courier_successor_variables
                                   for (arc1, arc2) in courier_transitions):
                                model.cbLazy(quicksum(self._courier_successor_variables[courier2, arc1, arc2]
                                                      for (arc1, arc2) in courier_transitions)
                                             <= len(courier_transitions) - 1)
                            # print(courier2, courier_transitions)
                        continue

                    # Iterate through the transitions, incrementing time
                    traversed_transitions: Set[Tuple[Arc, Arc]] = set()
                    while True:
                        # print(current_arc)
                        # Find the transition and the next arc
                        for (arc1, arc2) in courier_transitions:
                            if arc1 == current_arc:
                                # print(arc1, arc2)
                                # Found the transition, increment time
                                current_time += arc1.travel_time
                                traversed_transitions.add((arc1, arc2))
                                # Remove transition
                                if (arc1, arc2) in remaining_transitions:
                                    remaining_transitions.remove((arc1, arc2))

                                if current_time > arc2.latest_departure_time:
                                    # Infeasible set of arcs found, cut them
                                    # print('Infeasible time constraint')
                                    for courier2 in self._couriers:
                                        # Add a lazy constraint for all possible couriers
                                        if all((courier2, arc1, arc2) in self._courier_successor_variables
                                               for (arc1, arc2) in traversed_transitions):
                                            model.cbLazy(quicksum(
                                                self._courier_successor_variables[courier2, arc1, arc2]
                                                for (arc1, arc2) in traversed_transitions)
                                                         <= len(traversed_transitions) - 1)
                                            # print(courier2, traversed_transitions)
                                    current_time = arc2.earliest_departure_time
                                    traversed_transitions = set()
                                    current_arc = arc2

                                    # # Infeasible set of arcs found, trace backwards
                                    # TODO: Remove infinite loop
                                    # infeasible_transitions = trace_back_arcs(courier_transitions, arc2)
                                    # print('Infeasible time constraint')
                                    # for courier2 in self._couriers:
                                    #     # Add a lazy constraint for all possible couriers
                                    #     if all((courier2, arc1, arc2) in self._courier_successor_variables for
                                    #            (arc1, arc2) in infeasible_transitions):
                                    #         model.cbLazy(quicksum(
                                    #             self._courier_successor_variables[courier2, arc1, arc2] for (arc1, arc2)
                                    #             in infeasible_transitions) <= len(infeasible_transitions) - 1)
                                    #         print(courier2, infeasible_transitions)

                                else:
                                    # Change current arc to arc2, update current_time
                                    current_arc = arc2
                                    current_time = max(current_time, arc2.earliest_departure_time)

                                # Found the appropriate transition, move to find the next transition
                                break

                        # if type(arc2.arrival_location) == Group:
                        #     break

                        else:
                            # No transition found
                            assert False

                        if type(current_arc.arrival_location) == Group:
                            # Path followed to completion
                            # Remaining transitions form loops
                            if len(remaining_transitions) > 0:
                                # print('Cut off loops')
                                for courier2 in self._couriers:
                                    if all((courier2, arc1, arc2) in self._courier_successor_variables
                                           for (arc1, arc2) in remaining_transitions):
                                        model.cbLazy(quicksum(self._courier_successor_variables[courier2, arc1, arc2]
                                                              for (arc1, arc2) in remaining_transitions)
                                                     <= len(remaining_transitions) - 1)
                                        # print(courier2, remaining_transitions)
                            break

        super().optimize(callback)
        # print(f'Optimisation complete, t={math.ceil(time.time() - self._model_initiation)}')

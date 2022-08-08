#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 09:24:49 2022

@author: Tristan Hauber
"""

from gurobi import Model, quicksum, GRB
from classes import Group, Arc, Data, Order
from typing import List


class UntimedFragmentsMDRP(Model):
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
        """ ========== SETS ========== """
        self.couriers = group.couriers
        self.arcs = arcs
        self.orders = orders

        """ ========== VARIABLES ========== """
        # Payment to each courier
        self.courier_payments = {courier: self.addVar() for courier in self.couriers}
        # Time of departure for each arc
        self.timings = {arc: self.addVar() for arc in self.arcs}
        # If the order is delivered
        self.deliveries = {
            order: self.addVar(vtype=GRB.BINARY) for order in self.orders
        }
        # If the arc is serviced
        self.serviced = {arc: self.addVar(vtype=GRB.BINARY) for arc in self.arcs}
        # Courier assigned to arc
        self.assignments = {
            (courier, arc): self.addVar(vtype=GRB.BINARY)
            for courier in self.couriers
            for arc in self.arcs
        }
        # Arc2 follows arc1
        self.successors = {
            (arc1, arc2): self.addVar(vtype=GRB.BINARY)
            for arc1 in arcs
            for arc2 in arcs
            if arc2 != arc1
            and arc1.arrival_location == arc2.departure_location
            and arc1.earliest_departure_time + arc1.travel_time
            <= arc2.latest_departure_time
        }

        """ ========== OBJECTIVE ========== """
        self.setObjective(
            quicksum(self.courier_payments[courier] for courier in self.couriers)
            + quicksum(
                (1 - self.deliveries[order]) * cost_penalty for order in self.orders
            )
        )

        """ ========== CONSTRAINTS ========== """
        # Pay couriers for their time
        self.time_payments = {
            courier: self.addConstr(
                self.courier_payments[courier]
                >= (courier.off_time - courier.on_time) / 60 * Data.MIN_PAY_PER_HOUR
            )
            for courier in self.couriers
        }
        # Pay couriers for their deliveries
        self.delivery_payments = {
            courier: self.addConstr(
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
                order: self.addConstr(self.deliveries[order] == 1)
                for order in self.orders
            }
        # If delivered, an arc must service the order
        self.delivered = {
            order: self.addConstr(
                self.deliveries[order]
                == quicksum(
                    self.serviced[arc] for arc in self.arcs if order in arc.order_list
                )
            )
            for order in self.orders
        }
        # All activated arcs are assigned to a courier
        self.arcs_covered = {
            arc: self.addConstr(
                self.serviced[arc]
                == quicksum(
                    self.assignments[courier, arc]
                    for courier in self.couriers
                    if (courier, arc) in self.assignments
                )
            )
            for arc in self.arcs
        }
        # All activated non-exit arcs must have a successor
        predecessors = set()
        for (arc1, _) in self.successors:
            predecessors.add(arc1)
        self.has_successor = {
            arc1: self.addConstr(
                self.serviced[arc1]
                == quicksum(
                    self.successors[arc1, arc2]
                    for arc2 in arcs
                    if (arc1, arc2) in self.successors
                )
            )
            for arc1 in predecessors
        }
        # All successors must start late enough for the previous to finish
        self.succ_timings = {
            (arc1, arc2): self.addConstr(
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
            (arc1, arc2, courier): self.addConstr(
                self.assignments[courier, arc1] + self.successors[arc1, arc2] - 1
                <= self.assignments[courier, arc2]
            )
            for courier in self.couriers
            for (arc1, arc2) in self.successors
        }
        # If an arc is serviced, it must have a courier servicing it
        self.serv_requ_cour = {
            arc: self.addConstr(
                self.serviced[arc]
                == quicksum(self.assignments[courier, arc] for courier in self.couriers)
            )
            for arc in self.arcs
        }
        # An arc must be serviced after it is ready
        self.begin_on_time = {
            arc: self.addConstr(self.timings[arc] >= arc.earliest_departure_time)
            for arc in self.arcs
        }
        # An arc must be serviced before it is too late
        self.begin_in_time = {
            arc: self.addConstr(self.timings[arc] <= arc.latest_departure_time)
            for arc in self.arcs
        }

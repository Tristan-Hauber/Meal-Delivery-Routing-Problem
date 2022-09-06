from gurobipy import Model, GRB, quicksum
from classes import Data, Fragment, Courier, Order, Node
from typing import Set

class TimedFragmentsModel (Model):

    def __init__(self, fragments: Set[Fragment], couriers: Set[Courier], orders: Set[Order], nodes: Set[Node]) -> None:
        super().__init__()
        """ ========== Variables ========== """
        self._activations = {fragment: self.addVar(lb=0, ub=1, obj=0, vtype=GRB.BINARY, name=f'{fragment} activation',
                                                   column=None) for fragment in fragments}
        self._payments = {courier: self.addVar(lb=0, ub=1, obj=0, vtype=GRB.CONTINUOUS, name=f'{courier} payment',
                                               column=None) for courier in couriers}
        """ ========== Objective ========== """
        self.setObjective(quicksum(self._payments[courier] for courier in couriers))
        """ ========== Constraints ========== """
        self._delivery_payments = {courier: self.addConstr(self._payments[courier] >= quicksum(Data.PAY_PER_DELIVERY * len(fragment.order_list) * self._activations[fragment] for fragment in fragments if fragment.courier == courier), name=f'{courier} delivery payment') for courier in couriers}
        self._time_payments = {courier: self.addConstr(self._payments[courier] >= Data.MIN_PAY_PER_HOUR * (courier.off_time - courier.on_time) / 60, name=f'{courier} hourly payment') for courier in couriers}
        self._deliveries = {order: self.addConstr(quicksum(self._activations[fragment] for fragment in fragments if order in fragment.order_list) == 1, name=f'{order} delivery status') for order in orders}
        self._node_flow = {node: self.addConstr(quicksum(self._activations[fragment] for fragment in fragments if fragment.departure_node == node) == quicksum(self._activations[fragment] for fragment in fragments if fragment.arrival_node == node), name=f'{node} flow constraint') for node in nodes if type(node.location) != Courier}
        self._no_duplicate_couriers = {courier: self.addConstr(quicksum(self._activations[fragment] for fragment in fragments if type(fragment.departure_location) == Courier) <= 1, name=f'{courier} through network once') for courier in couriers}
        self._all_couriers_end = {courier: self.addConstr(quicksum(self._activations[fragment] for fragment in fragments if type(fragment.departure_location) == Courier) == quicksum(self._activations[fragment] for fragment in fragments if type(fragment.arrival_location) == Courier), name=f'{courier} must finish if start') for courier in couriers}
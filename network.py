from classes import Group, Order, Fragment
from typing import FrozenSet, Set

class SubNetwork:
    def __init__(self, group: Group, orders: FrozenSet[Order], objective: float, cumulative_cost_penalty: float, fragments: Set[Fragment]):
        self._group = group
        self._orders = orders
        self._best_objective_value = objective
        self._undelivered_orders_cost = cumulative_cost_penalty
        self._best_fragments = fragments

    def get_best_objective(self) -> float:
        return self._best_objective_value

    def get_best_fragments(self) -> Set[Fragment]:
        return self._best_fragments

    def get_group(self) -> Group:
        return self._group

    def get_orders(self) -> FrozenSet[Order]:
        return self._orders

    def get_undelivered_orders_cost(self) -> float:
        return self._undelivered_orders_cost
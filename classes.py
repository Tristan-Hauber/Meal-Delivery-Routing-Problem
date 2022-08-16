"""
A module with all classes needed to build and solve the MDRP.

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

Created on Mon Mar 14 16:43:08 2022

@author: Tristan
"""

from __future__ import annotations

import math

from typing import List, Set, Dict, Tuple


class Data:
    """Stores useful constants for the model.

    TRAVEL_SPEED: int
        The speed that a courier travels in metres per minute.
    PICKUP_SERVICE_TIME: int
        The time required to collect an order from a restaurant in minutes.
    DROPOFF_SERVICE_TIME: int
        The time required to drop off an order at a delivery locaation in minutes.
    TARGET_CLICK_TO_DOOR: int
        The targeted time between order placement and order delivery in minutes.
    MAX_CLICK_TO_DOOR: int
        The maximum time between order placement and order delivery in minutes.
    PAY_PER_DELIVERY: int
        The renumeration given to a courier for delivering an order in dollars.
    MIN_PAY_PER_HOUR: int
        The guaranteed minimum renumeration given to a courier for one hours worth of work in dollars per hour.
    """

    TRAVEL_SPEED = None
    PICKUP_SERVICE_TIME = None
    DROPOFF_SERVICE_TIME = None
    TARGET_CLICK_TO_DOOR = None
    MAX_CLICK_TO_DOOR = None
    PAY_PER_DELIVERY = None
    MIN_PAY_PER_HOUR = None


class Location:
    """A class containing data and functions to do with locations and distances.

    Properties
    ----------
    locations : List[Location]
        A list containing all created locations.
    distances : Dict[Tuple[Location, Location], float]
        A dictionary containing distances between pairs of locations.
    travel_times : Dict[Tuple[Location, Location], int]
        A dictionary containing travel times between pairs of locations.
    """

    locations = list()

    distances = dict()
    travel_times = dict()

    def __init__(self, x: int, y: int):
        """
        Create a new location with the given coordinates.

        Parameters
        ----------
        x : int
            The x grid coordinate of the location.
        y : int
            The y grid coordinate of the location.

        Returns
        -------
        None.

        """
        self.x = x
        self.y = y
        Location.locations.append(self)

    def get_distance_to(self, location: Location) -> float:
        """
        Return the distance to a given location.

        Parameters
        ----------
        location : Location
            The distination location.

        Returns
        -------
        float
            The distance to the destination location.

        """
        if (self, location) not in Location.distances:
            if (location, self) in Location.distances:
                Location.distances[(self, location)] = Location.distances[
                    (location, self)
                ]
            else:
                Location.distances[(self, location)] = math.sqrt(
                    (self.x - location.x) ** 2 + (self.y - location.y) ** 2
                )
        return Location.distances[(self, location)]

    def get_time_to(self, location: Location) -> int:
        """
        Return the time to a given location.

        Parameters
        ----------
        location : Location
            The destination location.

        Returns
        -------
        int
            The time required to travel to the destination location.

        """
        if (self, location) not in Location.travel_times:
            if (location, self) in Location.travel_times:
                Location.travel_times[(self, location)] = Location.travel_times[
                    (location, self)
                ]
            else:
                Location.travel_times[(self, location)] = math.ceil(
                    self.get_distance_to(location) / Data.TRAVEL_SPEED
                )
        return Location.travel_times[(self, location)]

    def get_location(self) -> Tuple[int, int]:
        """Return the location of the given object."""
        return (self.x, self.y)


class Courier(Location):
    """
    A class all about couriers.

    Class variables:
        couriers: List[Courier]
            A list containing all created couriers
        courier_by_code: Dict[str, Courier]
            A dictionary containing all couriers indexed by their unique code

    Instance variables:
        code: str
            The unique string referring to the courier
        x: int
            The x coordinate of the courier's starting location
        y: int
            The y coordinate of the courier's starting location
        on_time: int
            The time at which the courier's shift begins
        off_time: int
            The time at which the courier's shift ends
    """

    couriers = list()
    courier_by_code = dict()  # {str: Courier}

    def __init__(self, code: str, x: int, y: int, on_time: int, off_time: int):
        """
        Create a new Courier instance.

        Parameters
        ----------
        code : str
            The unique string referring to the courier.
        x : int
            The x coordinate of the courier's starting location.
        y : int
            The y coordinate of the courier's starting location.
        on_time : int
            The time that the courier's shift begins.
        off_time : int
            The time beyond which the courier may no longer pick up orders.

        Returns
        -------
        None.

        """
        super().__init__(x, y)
        self.code = code
        self.on_time = on_time
        self.off_time = off_time
        Courier.couriers.append(self)
        Courier.courier_by_code[self.code] = self

    def get_earliest_arrival_at(self, restaurant: Restaurant) -> int:
        """Return the earliest time a courier can arrive at a location."""
        return (
            self.on_time + self.get_time_to(restaurant) + Data.PICKUP_SERVICE_TIME / 2
        )

    def get_group(self) -> Group:
        """Return the group that the courier is part of."""
        for group in Group.groups:
            if self in group.couriers:
                return group
        # TODO: throw error (no group found)

    def __str__(self) -> str:
        """Return str(self)."""
        return self.code

    def __repr__(self) -> str:
        """Return repr(self)."""
        return self.__str__()


class Group(Location):
    """
    A class containing variables and functions to do with courier groups.

    Class variables:
        off_time: int
            The off_time of all couriers in the group
        couriers: [Courier]
            The list of all couriers in the group
        code: str
            The unique id of the group

    Instance variables:
        groups: [Group]
            A list of all created groups
        group_by_code: {str: Group}
            A dictionary containing all created groups, indexed by their id

    Class methods:
        group_couriers(Bool) -> None
            Group all of the couriers into courier groups
    """

    groups = list()

    def __init__(self, shift_end_time: int):
        """
        Create a new instance of Group.

        Parameters
        ----------
        shift_end_time : int
            The end time of all couriers in the group.

        Returns
        -------
        None.

        """
        super().__init__(0, 0)
        self.off_time = shift_end_time
        self.couriers = list()
        self.code = f"g{len(Group.groups) + 1}"
        self._arcs = None
        Group.groups.append(self)

    def get_earliest_arrival_at(self, restaurant: Restaurant) -> int:
        """Return the earliest arrival time at the given restaurant."""
        return min(
            courier.get_earliest_arrival_at(restaurant) for courier in self.couriers
        )

    def get_total_on_time(self) -> int:
        """Return the cumulative minutes worked by all couriers in the group."""
        return sum(courier.off_time - courier.on_time for courier in self.couriers)

    def group_couriers(group_by_off_time=True) -> None:
        """
        Group couriers into courier groups.

        Parameters
        ----------
        group_by_off_time : Bool
            True if all couriers that share an off_time are grouped, False otherwise.

        Returns
        -------
        None.

        """
        for courier in Courier.couriers:
            if group_by_off_time:
                found_group = False
                for group in Group.groups:
                    if group.off_time == courier.off_time:
                        group.couriers.append(courier)
                        found_group = True
                        break
                if not found_group:
                    group = Group(courier.off_time)
                    group.couriers.append(courier)
            else:  # Each courier is in its own group
                group = Group(courier.off_time)
                group.couriers.append(courier)

    def __str__(self) -> str:
        """Return str(self)."""
        return f"{self.code}"

    def __repr__(self) -> str:
        """Return repr(self)."""
        return self.__str__()

    def get_arcs(self) -> List[Arc]:
        """
        Retrieve a list of all arcs serviced by the group.

        Returns
        -------
        List[Arc]
            A list of arcs serviced by the group.

        """
        if self._arcs is None:
            self._arcs = list()
            for arc in Arc.arcs:
                if arc.group == self:
                    if arc.departure_location != arc.arrival_location or arc.order_list != []:
                        self._arcs.append(arc)
        return self._arcs


class Restaurant(Location):
    """
    A class containing variables and functions to do with restaurants.

    Class variables:
        restaurants: [Restaurant]
            A list containing all restaurants created
        restaurant_by_code: {str: Restaurant}
            A dictionary containing all restaurants indexed by their code

    Instance Variables:
        code: str
            The unique string referring to the restaurant
        x: int
            The x coordinate of the restaurant's location
        y: int
            The y coordinate of the restaurant's location
        orders: {Order}
            The set of orders placed at the restaurant

    Instance methods:
        get_orders() -> {Order}
            Return the set of orders placed at the restaurant
        get_latest_departure_time() -> int
            Return the latest latest_departure_time of orders at the restaurant
    """

    restaurants = list()
    restaurant_by_code = dict()

    def __init__(self, code: str, x: int, y: int):
        """
        Create a new instance of Restaurant.

        Parameters
        ----------
        code : str
            The unique string referring to the restaurant.
        x : int
            The x coordinate of the restaurant's location.
        y : int
            The y coordinate of the restaurant's location.

        Returns
        -------
        None.

        """
        super().__init__(x, y)
        self.code = code
        self.orders = set()

        Restaurant.restaurants.append(self)
        Restaurant.restaurant_by_code[code] = self

    def get_latest_departure_time(self) -> int:
        """Return the latest departure time of any order placed at this restaurant."""
        return max(order.latest_departure_time for order in self.orders)

    def __str__(self) -> str:
        """Return str(self)."""
        return self.code

    def __repr__(self) -> str:
        """Return repr(self)."""
        return self.__str__()


class Order(Location):
    """
    A class containing functions to do with Orders.

    Class variables:
        orders: [Order]
            A list containing all created orders
        order_by_code: {str: Order}
            A dictionary containing all orders indexed by their unique code

    Instance variables:
        code: str
            The unique code referring to the order
        x: int
            The x location of the order's dropoff location
        y: int
            The y location of the order's dropoff location
        placement_time: int
            The time at which the order was placed
        departure_restaurant: Restaurant
            The restaurant at which the order was placed
        ready_time: int
            The time after which the order can be collected
        latest_departure_time: int
            The time before which the order must be collected
    """

    orders = list()
    order_by_code = dict()  # {str: Order}

    def __init__(
        self,
        code: str,
        x: int,
        y: int,
        placement_time: int,
        restaurant: str,
        ready_time: int,
    ):
        """
        Create a new Order instance.

        Parameters
        ----------
        code : str
            The unique string referring to the order.
        x : int
            The x coordinate of the order's delivery location.
        y : int
            The y coordinate of the order's delivery location.
        placement_time : int
            The time at which the order was placed.
        restaurant : str
            The code for the restaurant at which the order was placed.
        ready_time : int
            The time at which the order can be picked up to be delivered.

        Returns
        -------
        None.

        """
        super().__init__(x, y)
        self.code = code
        self.placement_time = placement_time
        self.departure_restaurant = Restaurant.restaurant_by_code[restaurant]
        self.earliest_departure_time = ready_time
        Order.orders.append(self)
        Order.order_by_code[code] = self
        self.departure_restaurant.orders.add(self)
        self.latest_departure_time = (
            self.placement_time
            + Data.MAX_CLICK_TO_DOOR
            - self.departure_restaurant.get_time_to(self)
            - (Data.PICKUP_SERVICE_TIME + Data.DROPOFF_SERVICE_TIME) / 2
        )

    def __str__(self) -> str:
        """Return str(self)."""
        return self.code

    def __repr__(self) -> str:
        """Return repr(self)."""
        return self.__str__()

    def group_orders_by_restaurant(
        orders: List[Order],
    ) -> Dict[Restaurant, List[Order]]:
        """
        Group a list of orders by their departure restaurants.

        Parameters
        ----------
        orders : List[Order]
            The list of orders to group.

        Returns
        -------
        Dict[Restaurant, List[Order]]
            A dictionary of lists of orders, with the key being the common departure restaurant of the orders.

        """
        restaurants = dict()
        for order in orders:
            if order.departure_restaurant in restaurants:
                restaurants[order.departure_restaurant].append(order)
            else:
                restaurants[order.departure_restaurant] = [order]
        return restaurants


class Sequence:
    """
    A class containing variables and functions to do with order sequences.

    Class variables:
        sequences: [Sequence]
            A list containing all created sequences
        sequences_by_orders_and_last_order: {(frozenset(Order), Order): [Sequence]}
            A dictionary containing lists of sequences indexed by their order set and final order
        empty_sequences: [Sequence]
            A list containing all sequences with an empty order_list

    Instance variables:
        order_list: [Order]
            The list of orders serviced by this sequence
        ready_time: int
            The latest ready_time of any of the orders in order_list
        departure_location: Location
            The location from which the courier departs to deliver the orders
        total_travel_time: int
            The time taken to leave departure_location and deliver all orders in order_list
        latest_departure_time: int
            The latest time a courier can leave to deliver all orders in order_list on time

    Class methods:
        dominate([Sequence]) -> None
            Dominate the last sequence in the list against all other sequences in the list

    Instance methods:
        add_order(Order, dominate=True) -> Sequence
            Create a new sequence from the given sequence, adding Order to order_list
    """

    sequences = set()
    sequences_by_orders_and_last_order = dict()

    def __init__(self, order_list: List[Order], departure_location: Location):
        """
        Create a new instance of Sequence.

        The sequence does not necessarily have to contain orders; but if it does, then they must be placed at the restaurant given in departure_location. departure_location can be a restaurant or a courier.

        Parameters
        ----------
        order_list : List[Order]
            The list of orders, in the order to be delivered.
        departure_location : Location
            The location from which the courier departs.

        Returns
        -------
        None.

        """
        # TODO: if type(departure_location) != Courier and type(departure_location) != Restaurant: throw error
        self.order_list = order_list
        self.departure_location = departure_location
        # TODO: if len(self.order_list) > 0: if self.order_list[0].get_departure_restaurant() != self.departure_location: throw error

        # Set earliest_departure_time
        if len(self.order_list) > 0:
            # earliest_departure_time is the latest ready time of the orders
            self.earliest_departure_time = max(
                order.earliest_departure_time for order in self.order_list
            )
        else:
            # no orders, check location for earliest_departure_time
            if type(self.departure_location) == Restaurant:
                # earliest_departure_time is the earliest ready time of the orders
                self.earliest_departure_time = min(
                    order.earliest_departure_time
                    for order in self.departure_location.orders
                )
            else:
                # earliest_departure_time is the on time of the courier
                self.earliest_departure_time = self.departure_location.on_time

        # Set travel_time and latest_departure_time
        if len(self.order_list) > 0:
            # Iterate through the orders, calculating travel_time and latest_departure_time
            self.travel_time = (
                self.departure_location.get_time_to(self.order_list[0])
                + (Data.PICKUP_SERVICE_TIME + Data.DROPOFF_SERVICE_TIME) / 2
            )
            self.latest_departure_time = (
                self.order_list[0].placement_time
                + Data.MAX_CLICK_TO_DOOR
                - self.travel_time
            )

            for i in range(len(self.order_list) - 1):
                # Iterate through all orders except the last one, adding the journey to the total
                self.travel_time += (
                    self.order_list[i].get_time_to(self.order_list[i + 1])
                    + Data.DROPOFF_SERVICE_TIME
                )
                self.latest_departure_time = min(
                    self.latest_departure_time,
                    self.order_list[i + 1].placement_time
                    + Data.MAX_CLICK_TO_DOOR
                    - self.travel_time,
                )

        else:
            # No orders in sequence. travel_time is zero
            self.travel_time = 0

            if type(self.departure_location) == Restaurant:
                # latest_departure_time is the latest latest_departure_time of any order at the restaurant
                self.latest_departure_time = max(
                    order.latest_departure_time
                    for order in self.departure_location.orders
                )

            else:
                # latest_departure_time is the courier's off time
                self.latest_departure_time = self.departure_location.off_time

        Sequence.sequences.add(self)
        # Dominate the sequence
        if len(self.order_list) > 2:
            dominate_key = (frozenset(self.order_list), self.order_list[-1])
            if dominate_key not in Sequence.sequences_by_orders_and_last_order:
                Sequence.sequences_by_orders_and_last_order[dominate_key] = [self]
            else:
                dominate_list = Sequence.sequences_by_orders_and_last_order[
                    dominate_key
                ]
                dominate_list.append(self)
                dominated_sequences = Sequence.dominate(dominate_list)
                for sequence in dominated_sequences:
                    dominate_list.remove(sequence)
                    Sequence.sequences.remove(sequence)

    def add_order(self, order: Order) -> Sequence:
        """
        Create a new sequence by adding an order to the current sequence.

        Parameters
        ----------
        order : Order
            The order to add to the sequence.

        Returns
        -------
        Sequence
            The sequence with the order added.

        """
        # TODO: if order in self.order_list: throw error
        new_list = list(self.order_list)
        new_list.append(order)
        return Sequence(new_list, self.departure_location)

    def dominate(sequences: List[Sequence]) -> Set[Sequence]:
        """
        Dominate the given sequences.

        Given two sequences s1 and s2, s1 dominates s2 if s1 has a better or
        equal travel time and latest leaving time. s2 dominates s1 if s1
        doesn't dominate s2 and s2 has a better or equal travel time and
        latest leaving time.

        Assume that the sequence to be dominated is the last one of the list.
        This function compares all sequences in the list with the last one,
        keeping track of all dominated sequences. It then removes the
        dominated sequences from all sequence lists.

        Parameters
        ----------
        sequences : List[Sequence]
            A list of sequences to be dominated.

        Returns
        -------
        dominated_sequences : Set[Sequence]

        """
        dominated_sequences = set()
        sequence_one = sequences[-1]
        for sequence_two in sequences[:-1]:
            if (
                sequence_one.travel_time <= sequence_two.travel_time
                and sequence_one.latest_departure_time
                >= sequence_two.latest_departure_time
            ):
                dominated_sequences.add(sequence_two)
            elif (
                sequence_two.travel_time <= sequence_one.travel_time
                and sequence_two.latest_departure_time
                >= sequence_one.latest_departure_time
            ):
                dominated_sequences.add(sequence_one)
        return dominated_sequences

    def __str__(self) -> str:
        """Return str(self)."""
        return f"({self.departure_location}, {self.order_list})"

    def __repr__(self) -> str:
        """Return repr(self)."""
        return self.__str__()


class Arc:
    """
    A class containing variables and functions to do with untimed path fragments.

    Class variables:
        arcs: [Arc]
            A list of all created arcs
        arcs_by_group_and_order_set_and_arrival_location: {(Group, frozenset(Order), Location): [Arc]}
            A dictionary containing lists of arcs indexed by the servicing courier, serviced orders, and next location
        arcs_by_group_and_arrival_location: {(Group, Location): [Arc]}
            A dictionary containing lists of arcs indexed the the servicing group and arrival location

    Instance variables:
        group: Group
            The group servicing the arc
        sequence: Sequence
            The order sequence serviced by the arc
        arrival_location: Location
            The location the courier goes to after servicing the arc
        order_list: [Order]
            The orders in sequence serviced by the arc
        departure_location: Location
            The location the courier departs from to service the arc
        total_travel_time: int
            The time taken to depart from departure_location, service sequence, and arrive at arrival_location
        earliest_departure_time: int
            The earliest courier can leave departure_location given the time taken to get their and sequence.ready_time
        earliest_arrival_time: int
            The earliest courier can arrive at arrival_location
        latest_arrival_time: int
            The latest courier can arrive at arrival_location to pick up a new order before their shift ends
        latest_departure_time: int
            The latest courier can leave departure_location given sequence.latest_departure_time and latest_arrival_time

    Class methods:
        dominate([Arc]) -> {Arc}
            Dominate the given list of arcs and return all dominated arcs.
    """

    arcs = set()
    arcs_by_group_and_order_set_and_arrival_location = dict()
    arcs_by_group_and_arrival_location = dict()
    arcs_by_group_and_departure_location = dict()
    arc_by_components = dict()

    def __init__(self, group: Group, sequence: Sequence, arrival_location: Location):
        # Given attributes
        self.group = group
        self.sequence = sequence
        self.arrival_location = arrival_location

        # Computed attributes
        self.order_list = self.sequence.order_list
        self.departure_location = self.sequence.departure_location

        # earliest_departure_time
        if type(self.departure_location) == Restaurant:
            self.earliest_departure_time = max(
                self.group.get_earliest_arrival_at(self.departure_location),
                self.sequence.earliest_departure_time,
            )
        else:  # type(self.departure_location) == Courier
            assert type(self.departure_location) == Courier, "incorrect location type"
            self.earliest_departure_time = self.departure_location.on_time

        # travel_time
        if len(self.sequence.order_list) > 0:
            if type(self.arrival_location) == Restaurant:  # Main arc
                last_order = self.sequence.order_list[-1]
                self.travel_time = (
                    self.sequence.travel_time
                    + last_order.get_time_to(self.arrival_location)
                    + (Data.PICKUP_SERVICE_TIME + Data.DROPOFF_SERVICE_TIME) / 2
                )
            else:  # Departure arc
                assert type(self.arrival_location) == Group, "incorrect location type"
                self.travel_time = self.sequence.travel_time
        else:
            if type(self.departure_location) == Courier:  # Entry arc
                assert (
                    type(self.arrival_location) == Restaurant
                ), "incorrect location type"
                self.travel_time = (
                    self.departure_location.get_time_to(self.arrival_location)
                    + Data.PICKUP_SERVICE_TIME / 2
                )
            else:  # Waiting arc
                assert (
                    type(self.departure_location) == Restaurant
                ), "incorrect location type"
                self.travel_time = 1

        # earliest_arrival_time
        self.earliest_arrival_time = self.earliest_departure_time + self.travel_time

        # latest_arrival_time
        if type(self.arrival_location) == Restaurant:
            deliverable_orders = set(
                order
                for order in self.arrival_location.orders
                if order.earliest_departure_time <= self.group.off_time
                if order not in self.sequence.order_list
            )
            assert len(deliverable_orders) > 0
            self.latest_arrival_time = min(
                self.group.off_time,
                max(order.latest_departure_time for order in deliverable_orders),
            )
        else:
            assert type(self.arrival_location) == Group
            self.latest_arrival_time = (
                min(self.group.off_time, self.sequence.latest_departure_time)
                + self.travel_time
            )

        # latest_departure_time
        self.latest_departure_time = min(
            self.group.off_time,
            self.sequence.latest_departure_time,
            self.latest_arrival_time - self.travel_time,
        )
        # update latest_arrival_time
        self.latest_arrival_time = min(
            self.latest_arrival_time, self.latest_departure_time + self.travel_time
        )

        # Dominate
        if self.earliest_departure_time <= self.latest_departure_time:
            Arc.arcs.add(self)
            key = (self.group, self.arrival_location)
            dominate_key = (
                self.group,
                frozenset(self.sequence.order_list),
                self.arrival_location,
            )
            if key not in Arc.arcs_by_group_and_arrival_location:
                Arc.arcs_by_group_and_arrival_location[key] = []
            if dominate_key not in Arc.arcs_by_group_and_order_set_and_arrival_location:
                Arc.arcs_by_group_and_order_set_and_arrival_location[dominate_key] = []
            if (
                self.group,
                self.departure_location,
            ) not in Arc.arcs_by_group_and_departure_location:
                Arc.arcs_by_group_and_departure_location[
                    (self.group, self.departure_location)
                ] = list()
            Arc.arcs_by_group_and_arrival_location[key].append(self)
            Arc.arcs_by_group_and_order_set_and_arrival_location[dominate_key].append(
                self
            )
            Arc.arcs_by_group_and_departure_location[
                self.group, self.departure_location
            ].append(self)
            if len(self.sequence.order_list) >= 2:
                dominate_list = Arc.arcs_by_group_and_order_set_and_arrival_location[
                    dominate_key
                ]
                if len(dominate_list) >= 2:
                    dominated_arcs = Arc.dominate(dominate_list)
                    for arc in dominated_arcs:
                        dominate_list.remove(arc)
                        Arc.arcs.remove(arc)
                        Arc.arcs_by_group_and_arrival_location[
                            arc.group, arc.arrival_location
                        ].remove(arc)
                        Arc.arcs_by_group_and_departure_location[
                            arc.group, arc.departure_location
                        ].remove(arc)
        if self in Arc.arcs:
            Arc.arc_by_components[
                (self.group, self.sequence, self.arrival_location)
            ] = self

    def dominate(arcs: List[Arc]) -> Set[Arc]:
        """
        Dominate a list of untimed path fragments.

        The function takes as input a list of untimed path fragments. The last fragment in this
        list is to be dominated against all other fragments in the list. We assume that all
        other fragments in the list share an order list, and departure and arrival location.
        Then one fragment dominates another if it has a later (or equal) departure time and
        a shorter (or equal) travel time.

        Parameters
        ----------
        arcs : List[Arc]
            The list of untimed path fragments to be dominated.

        Returns
        -------
        Set[Arc]
            The set of untimed path fragments after domination.

        """
        arc1 = arcs[-1]
        dominated_arcs = set()
        for arc2 in arcs[:-1]:
            if (
                arc1.latest_departure_time >= arc2.latest_departure_time
                and arc1.travel_time <= arc2.travel_time
            ):
                dominated_arcs.add(arc2)
            elif (
                arc2.latest_departure_time >= arc1.latest_departure_time
                and arc2.travel_time <= arc1.travel_time
            ):
                dominated_arcs.add(arc1)
        return dominated_arcs

    def get_pred(self) -> Set[Arc]:
        """
        Return a set of predecessors to the given arc.

        Parameters
        ----------
        arc : Arc
            The arc to get the set of predecessors for.

        Returns
        -------
        predecessors : set(Arc)
            All arcs that are predecessors to the given arc.

        """
        if type(self.departure_location) == Courier:
            return set()
        assert type(self.departure_location) == Restaurant
        departure_location = self.departure_location
        latest_departure_time = self.latest_departure_time
        group = self.group
        return set(
            arc
            for arc in Arc.arcs_by_group_and_arrival_location[
                (group, departure_location)
            ]
            if arc.earliest_arrival_time <= latest_departure_time
            if set(arc.sequence.order_list).intersection(set(self.sequence.order_list))
            == set()
            if (
                arc.departure_location != arc.arrival_location
                or arc.sequence.order_list != []
            )
        )

    def get_succ(self) -> Set[Arc]:
        """Find and return all possible successors of the given untimed path fragment."""
        if type(self.arrival_location) == Group:
            return set()
        assert type(self.arrival_location) == Restaurant
        restaurant: Restaurant = self.arrival_location
        earliest_arrival_time: int = self.earliest_arrival_time
        group: Group = self.group
        succ_arcs: Set[Arc] = set()
        for arc in Arc.arcs_by_group_and_departure_location[group, restaurant]:
            if set(self.order_list).intersection(set(arc.order_list)) == set():
                if arc.latest_departure_time >= earliest_arrival_time:
                    if arc.order_list != []:
                        succ_arcs.add(arc)
        return succ_arcs

    def get_pred_to_arcs(arcs: set(Arc)) -> set(Arc):
        """
        Return a set of predecessors to a set of arcs.

        Parameters
        ----------
        arcs : set(Arc)
            The set of arcs to get the set of predecessors for.

        Returns
        -------
        predecessors : set(Arc)
            The set of predecessors to the given arcs.

        """
        predecessors = set()
        for arc in arcs:
            predecessors.union(arc.get_pred())
        return predecessors

    def get_arc_from_components(
        group: Group, sequence: Sequence, arrival_location: Location
    ) -> Arc:
        """Find and return the unique arc with the given components."""
        return Arc.arc_by_components[(group, sequence, arrival_location)]

    def __str__(self) -> str:
        """Return str(self)."""
        return f"({self.group}, {self.sequence}, {self.arrival_location})"

    def __repr__(self) -> str:
        """Return repr(self)."""
        return self.__str__()

    def get_arcs_for_orders(orders: List[Order], group: Group) -> List[Arc]:
        """
        Get all arcs that deliver orders from the given list, and are serviced by the given group.

        Parameters
        ----------
        orders : List[Order]
            The list of orders to retrieve arcs for.
        group : Group
            The courier group to retrieve arcs for.

        Returns
        -------
        List[Arc]
            The list of arcs containing the given orders, serviced by the given group.

        """
        order_set = set(orders)

        arrival_locations = list()
        for restaurant in Order.group_orders_by_restaurant(orders):
            arrival_locations.append(restaurant)
        arrival_locations.append(group)

        arcs = list()
        for key in Arc.arcs_by_group_and_order_set_and_arrival_location.keys():
            if (
                key[0] == group
                and key[1].issubset(order_set)
                and key[2] in arrival_locations
            ):
                for arc in Arc.arcs_by_group_and_order_set_and_arrival_location[key]:
                    arcs.append(arc)

        return arcs


class Node:
    """
    A class containing information about points in group-location-time space.

    Class variables:
        nodes: [Node]
            A list containing all created nodes.
        node_by_components: {(Group, Location, int): Node}
            A dictionary containing all nodes indexed by the group, location and time
        nodes_by_group_and_location: {(Group, Location): [Node]}
            A dictionary containing lists of nodes indexed by the courier and location

    Instance variables:
        courier: Courier
            The courier sub-network the node is part of
        location: Location
            The location of the node
        time: int
            The time location of the node
    """

    nodes = list()
    node_by_components = dict()
    nodes_by_group_and_location = dict()

    def __init__(self, group, location, time):
        """Create a new point in GLT space."""
        self.group = group
        self.location = location
        self.time = time
        Node.nodes.append(self)
        Node.node_by_components[(group, location, time)] = self
        if (group, location) not in Node.nodes_by_group_and_location:
            Node.nodes_by_group_and_location[(group, location)] = list()
        Node.nodes_by_group_and_location[(group, location)].append(self)

    def get_arrival_node(group: Group, restaurant: Restaurant, time: int) -> Node:
        """
        Return a node based on the given group, restaurant and time.

        Round the time down if there are deliverable orders, up otherwise.

        Parameters
        ----------
        group : Group
            The group arriving at the node.
        restaurant : Restaurant
            The arrival restaurant to find a node for.
        time : int
            The arrival time that needs to be rounded.

        Returns
        -------
        Node
            The arrival node, after rounding the time appropriately.

        """
        # Check to see if node at arrival time
        if (group, restaurant, time) in Node.node_by_components:
            return Node.node_by_components[(group, restaurant, time)]

        # No node at arrival time
        # Check if active order at arrival time
        # If active order, round down. Otherwise, round up.
        active_orders = list(
            order
            for order in restaurant.orders
            if order.earliest_departure_time <= time
            if order.latest_departure_time >= time
        )
        if len(active_orders) > 0:
            arrival_time = max(
                node.time
                for node in Node.nodes_by_group_and_location[(group, restaurant)]
                if node.time < time
            )
        else:
            arrival_time = min(
                node.time
                for node in Node.nodes_by_group_and_location[(group, restaurant)]
                if node.time > time
            )

        # Return the calculated node
        return Node.node_by_components[(group, restaurant, arrival_time)]

    def __str__(self):
        """Return str(self)."""
        return f"({self.group}, {self.location}, {self.time})"

    def __repr__(self):
        """Return repr(self)."""
        return self.__str__()


class Fragment:
    """
    A class containing variables and functions to do with timed path fragments.

    Class variables:
        fragments: [Fragment]
            A list containing all created path fragments.
        fragments_by_departure_node: {Node: [Fragment]}
            A dictionary containing lists of fragments indexed by the departure node
        fragments_by_arrival_node: {Node: [Fragment]}
            A dictionary containing lists of frgaments indexed by the arrival node
        fragments_by_sequence: {Sequence: [Fragment]}
            A dictionary containing lists of fragments indexed by the serviced sequence
        fragments_by_arc: {Arc: [Fragment]}
            A dictionary containing lists of fragments indexed by the serviced arc
        fragments_by_group: {Group: [Fragment]}
            A dictionary containing lists of fragments indexed by the servicing group
        fragments_by_order: {Order: [Fragment]}
            A dictionary containing lists of fragments indexed by a contained order
        departure_fragments_by_group: {Group: [Fragment]}
            A dictionary containing lists of fragments departing from depots indexed by the servicing courier

    Instance variables:
        departure_node: Node
            The node from which the fragment departs
        arc: Arc
            The arc serviced by the fragment
        group: Group
            The group servicing the fragment
        departure_location: Location
            The location from which the fragment departs
        departure_time: int
            The time at which the fragment departs
        sequence: Sequence
            The sequence serviced by the fragment
        arrival_location: Location
            The location at which the fragment arrives
        order_list: [Order]
            The orders serviced by the fragment
        arrival_time: int
            The time at which the fragment arrives
        arrival_node: Node
            The node at which the fragment arrives

    Class methods:
        dominate([Fragment]) -> {Fragment}
            Dominate the given list of fragments, and return a set of dominated fragments.
        sort() -> None
            Sort all fragments into useful dictionaries for later accessibility.
    """

    fragments = list()
    fragments_by_departure_node = dict()
    fragments_by_arrival_node = dict()
    fragments_by_sequence = dict()
    fragments_by_arc = dict()
    fragments_by_group = dict()
    fragments_by_order = dict()
    departure_fragments_by_courier = dict()

    def __init__(self, departure_node, arc):
        """
        Create a new instance of Fragment.

        Parameters
        ----------
        departure_node : Node
            The node from which the fragment departs.
        arc : Arc
            The arc serviced by the fragment.
        dominate : Bool, optional
            True if dominate the fragments, false otherwise. The default is True.

        Returns
        -------
        None.

        """
        self.departure_node = departure_node
        self.arc = arc
        self.group = arc.group
        self.departure_location = arc.departure_location
        self.departure_time = departure_node.time
        self.sequence = arc.sequence
        self.arrival_location = arc.arrival_location
        self.order_list = self.arc.sequence.order_list

        # Find the arrival node for the arc
        if type(self.arrival_location) == Group:
            self.arrival_node = Node.node_by_components[
                (self.group, self.arrival_location, self.group.off_time)
            ]
        elif type(self.arrival_location) == Restaurant:
            # If waiting arc, move to the next time
            # If not waiting arc, use calculated arrival node
            if (
                self.order_list == []
                and self.departure_location == self.arrival_location
            ):
                # Waiting path fragment, increment time
                arrival_time = min(
                    node.time
                    for node in Node.nodes_by_group_and_location[
                        (self.group, self.arrival_location)
                    ]
                    if node.time > self.departure_time
                )
                self.arrival_node = Node.node_by_components[
                    (self.group, self.arrival_location, arrival_time)
                ]
            else:
                # Use calculated arrival_node
                self.arrival_node = Node.get_arrival_node(
                    self.group,
                    self.arrival_location,
                    self.departure_time + self.arc.travel_time,
                )
        else:
            assert False

        # Set arrival time
        self.arrival_time = self.arrival_node.time

        # Try dominating against the last fragment added, rather than all arcs so far added.
        if len(Fragment.fragments) > 0:
            last_fragment = Fragment.fragments[-1]
            if (
                self.arc == last_fragment.arc
                and self.arrival_time == last_fragment.arrival_time
                and self.departure_time > last_fragment.departure_time
                and self.departure_location == last_fragment.departure_location
            ):
                Fragment.fragments[-1] = self
            else:
                Fragment.fragments.append(self)
        else:
            Fragment.fragments.append(self)

    def dominate(fragments):
        """
        Dominate the given list of fragments.

        The method will compare each fragment in the list with the last member,
        and return a set of all dominated fragments.

        Parameters
        ----------
        fragments : [Fragment]
            The fragments to dominate.

        Returns
        -------
        dominated_fragments : {Fragment}
            A set of all dominated fragments.

        """
        fragment_one = fragments[-1]
        dominated_fragments = set()
        for fragment_two in fragments[:-1]:
            if (
                fragment_one.departure_node.time >= fragment_two.departure_node.time
                and fragment_one.arrival_node.time <= fragment_two.arrival_node.time
            ):
                dominated_fragments.add(fragment_two)
            elif (
                fragment_two.departure_node.time >= fragment_one.departure_node.time
                and fragment_two.arrival_node.time <= fragment_one.arrival_node.time
            ):
                dominated_fragments.add(fragment_one)
        return dominated_fragments

    def sort():
        """
        Sort all fragments into useful dictionaries.

        Returns
        -------
        None.

        """
        for node in Node.nodes:
            Fragment.fragments_by_departure_node[node] = list()
            Fragment.fragments_by_arrival_node[node] = list()

        for sequence in Sequence.sequences:
            Fragment.fragments_by_sequence[sequence] = list()

        for arc in Arc.arcs:
            Fragment.fragments_by_arc[arc] = list()

        for group in Group.groups:
            Fragment.fragments_by_group[group] = list()

        for courier in Courier.couriers:
            Fragment.departure_fragments_by_courier[courier] = list()

        for order in Order.orders:
            Fragment.fragments_by_order[order] = list()

        for fragment in Fragment.fragments:
            Fragment.fragments_by_departure_node[fragment.departure_node].append(
                fragment
            )
            Fragment.fragments_by_arrival_node[fragment.arrival_node].append(fragment)

            Fragment.fragments_by_sequence[fragment.sequence].append(fragment)

            if fragment.arc not in Fragment.fragments_by_arc:
                Fragment.fragments_by_arc[fragment.arc] = list()
            Fragment.fragments_by_arc[fragment.arc].append(fragment)

            Fragment.fragments_by_group[fragment.group].append(fragment)

            for order in fragment.sequence.order_list:
                Fragment.fragments_by_order[order].append(fragment)

            if type(fragment.departure_location) == Courier:
                courier: Courier = fragment.departure_location
                Fragment.departure_fragments_by_courier[courier].append(fragment)

    def get_fragment_by_arc_and_departure_node(
        arc: Arc, departure_node: Node
    ) -> Fragment:
        """
        Return the unique fragment defined by the given arc and departure node.

        Parameters
        ----------
        arc : Arc
            The arc that the fragment is built on.
        departure_node : Node
            The node the arc departs from.

        Returns
        -------
        Fragment
            The fragment containing arc and departing from departure_node.

        """
        for fragment in Fragment.fragments_by_arc[arc]:
            if fragment.departure_node == departure_node:
                return fragment

    def __str__(self):
        """Return str(self)."""
        return f"({self.departure_node}, {self.arc}, {self.arrival_node})"

    def __repr__(self):
        """Return repr(self)."""
        return self.__str__()

    def get_timed_fragments_from_untimed_fragment_network(
        paths: List[UntimedFragmentPath],
    ) -> List[Fragment]:
        """
        Create a list of timed fragments from a list of untimed fragment paths.

        This function takes as input a list of untimed fragment paths representing a complete untimed fragment network, and converts each path to one consisting of timed fragments. It then collects all the timed fragments in each path and returns them as a completed list.

        Parameters
        ----------
        paths : List[UntimedFragmentPath]
            A list of untimed fragment paths.

        Returns
        -------
        List[Fragment]
            A list of timed fragments built on the untimed fragment paths.

        """
        timed_path_fragments = list()
        for path in paths:
            timed_fragment_path = path.convert_to_timed_fragment_path()
            for timed_fragment in timed_fragment_path.path:
                timed_path_fragments.append(timed_fragment)
        return timed_path_fragments

    def get_waiting_timed_fragment_from_node(node: Node) -> Fragment:
        """
        Find the (unique) waiting timed fragment from a given node.

        Parameters
        ----------
        node : Node
            The node to search for a waiting timed fragment from.

        Returns
        -------
        timed_fragment : Fragment
            The waiting timed fragment departing from the given node.

        """
        for timed_fragment in Fragment.fragments_by_departure_node[node]:
            if timed_fragment.order_list == list() and timed_fragment.arrival_location == timed_fragment.departure_location:
                return timed_fragment

    def get_fragments_from_orders(orders: Set[Order]) -> Set[Fragment]:
        """Find and return the set of fragments that deliver the given set of orders."""
        fragments = set()
        for order in orders:
            fragments.union(Fragment.fragments_by_order[order])
        return fragments


class UntimedFragmentPath:
    """
    A class representing a string of untimed fragments joined together to make a complete path through the network.

    Attributes
    ----------
    path : List[Arc]
        The untimed fragments that make up the path through the network.
    """

    def __init__(self, path: List[Arc]):
        self.path = path

    def get_untimed_fragments(self) -> List[Arc]:
        """
        Get the untimed fragments that make up the path as a list.

        Returns
        -------
        List[Arc]
            The list of untimed fragments that make up the path.

        Methods
        -------
        convert_to_timed_fragment_path(self) -> TimedFragmentPath
            Convert the path of untimed fragments into a path of timed fragments.

        """
        return self.path

    def convert_to_timed_fragment_path(self) -> TimedFragmentPath:
        """
        Convert the path from being made up of untimed fragments to one being made up of timed fragments.

        Returns
        -------
        TimedFragmentPath
            The timed fragment path equivalent to the untimed fragment path.

        """
        timed_fragment_path = list()
        current_time = self.path[0].earliest_departure_time
        for untimed_fragment in self.path:
            # Iterate through all the untimed fragments
            current_time = min(current_time, untimed_fragment.earliest_departure_time)
            best_timed_fragment = None
            for timed_fragment in Fragment.fragments_by_arc[untimed_fragment]:
                # Find the earliest timed fragment
                if timed_fragment.departure_time >= current_time and (best_timed_fragment is None or timed_fragment.departure_time < best_timed_fragment.departure_time):
                    best_timed_fragment = timed_fragment
            if best_timed_fragment.departure_time != current_time:
                # We need to add waiting arcs
                group = best_timed_fragment.group
                location = best_timed_fragment.departure_location
                while current_time < best_timed_fragment.departure_time:
                    # Choose the node with the largest time before the current time (round the time down). If no nodes, round the time up to the lowest time.
                    possible_nodes = list(node for node in Node.nodes_by_group_and_location[(group, location)] if node.time <= current_time)
                    possible_nodes.sort(key=lambda node: node.time)
                    if len(possible_nodes) == 0:
                        possible_nodes = list(node for node in Node.nodes_by_group_and_location[(group, location)] if node.time > current_time)
                        possible_nodes.sort(key=lambda node: node.time)
                        node = possible_nodes[0]
                    else:
                        node = possible_nodes[-1]
                    waiting_timed_fragment = Fragment.get_waiting_timed_fragment_from_node(node)
                    timed_fragment_path.append(waiting_timed_fragment)
                    current_time = waiting_timed_fragment.arrival_time
            # Add fragment to the list
            timed_fragment_path.append(best_timed_fragment)
            # Update time
            current_time = best_timed_fragment.arrival_time
        return TimedFragmentPath(timed_fragment_path)


class TimedFragmentPath:
    """
    A class representing a string of timed fragments joined together to make a complete path through the network.

    Attributes
    ----------
    path : List[Fragment]
        The list of timed fragments that make up the path.
    """

    def __init__(self, path: List[Fragment]):
        self.path = path

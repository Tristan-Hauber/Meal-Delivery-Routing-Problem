"""A short script to check if a given solution is valid."""
from classes import *
from typing import *


def check_solution(orders: Set[Order], itineraries: List[UntimedFragmentPath]) -> bool:
    """
    A short function to check if a given solution is valid.

    A valid solution is one that:
    - Delivers all orders it claims to deliver
    - All orders are delivered on time
    - All orders are picked up within the couriers' shifts
    - All orders are picked up from the placement restaurant
    - All orders are delivered at most once
    - All orders are picked up after their ready time
    """

    serviced_orders: Set[Order] = set()
    for itinerary in itineraries:
        arcs: List[Arc] = itinerary.path
        courier: Courier = arcs[0].departure_location
        time = courier.on_time
        location: Location = courier
        for arc in arcs[1:]:
            time += location.get_time_to(arc.departure_location) + Data.PICKUP_SERVICE_TIME / 2
            earliest_pickup = max(order.earliest_departure_time for order in arc.order_list)
            time = max(time, earliest_pickup)
            if time > courier.off_time:
                print(f'{courier} arrived at {arc.departure_location} after off time')
                return False
            location = arc.departure_location
            time += Data.PICKUP_SERVICE_TIME / 2
            for order in arc.order_list:
                if order in serviced_orders:
                    print(f'{order} delivered at least twice')
                    return False
                if order.departure_restaurant != arc.departure_location:
                    print(f'{order} picked up from the wrong location')
                    return False
                serviced_orders.add(order)
                time += location.get_time_to(order) + Data.DROPOFF_SERVICE_TIME / 2
                if time > order.placement_time + Data.MAX_CLICK_TO_DOOR:
                    print(f'{order} not delivered in time')
                    return False
                location = order
                time += Data.DROPOFF_SERVICE_TIME / 2

    if serviced_orders != orders:
        print(f'Did not deliver orders {orders.difference(serviced_orders)}')
        return False

    print(f'All orders serviced on time, and within courier shifts.')
    return True

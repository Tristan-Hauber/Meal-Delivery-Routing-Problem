#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module dedicated to reading and interpreting a model solution.

Created on Sat Apr 23 14:50:46 2022

@author: Tristan

Classes
-------
Assignment
    A class holding information about courier assignments.
Delivery
    A class holding information about order deliveries.
Dispatch
    A class holding information about individual courier dispatches.

Attributes
----------
assignment_information: [Assignment]
    A list of all courier assignments.
order_delivery_information: [Delivery]
    A list of all order deliveries.
courier_dispatch_information: [Dispatch]
    A list of all courier dispatches.

"""

from classes import Location, Data, Restaurant, Fragment, Order, Courier, Group


class Assignment:
    """A class holding information about courier assignments.

    A courier assignment can be summarised as an assignment time, or when the
    courier was given the assignment, a pickup time, i.e. the time the courier
    picks up the orders in the assignment, a courier ID, to identify the
    courier given the assignment, and a list of order IDs, to identify the
    orders in the assignment, and the order in which they are delivered.

    Attributes
    ----------
    assignment_time: int
        The time the assignment was given to the courier
    pickup_time: int
        The time the courier picked up the assignment.
    courier_ID: str
        The ID of the courier that was given the assignment.
    order_IDs: [str]
        A list of IDs of all orders in the assignment.

    """

    def __init__(self, fragment: Fragment) -> None:
        """
        Create a new courier assignment.

        Parameters
        ----------
        fragment : Fragment
            A path fragment from which to create the courier assignment.

            fragment must have at least one order, that is
            len(fragment.order_list) > 0.

        Returns
        -------
        None

        """
        self.assignment_time: int = fragment.departure_time
        self.pickup_time: int = fragment.departure_time
        self.courier_ID: str = fragment.group.couriers[0].code
        self.order_IDs: [str] = list(order.code for order in
                                     fragment.arc.sequence.order_list)

    def __str__(self) -> str:
        """
        Format the assignment into a string.

        The string returned always comes in the following format:
        assignment time, pickup time, Courier ID, Order ID, ..., Order ID

        Returns
        -------
        str
            A text description of the assignment.

        """
        string = f'{self.assignment_time} {self.pickup_time} {self.courier_ID}'
        for ID in self.order_IDs:
            string += f' {ID}'
        return string


class Delivery:
    """A class holding information about order deliveries.

    An order delivery is a set of information including:
        - Order ID
        - Placement time
        - Ready time
        - Pickup time
        - Delivery time
        - Courier ID

    Attributes
    ----------
    order_ID: str
        The string ID of the order.
    placement_time: int
        The time at which the order was placed.
    ready_time: int
        The time at which the order is ready for pickup.
    pickup_time: int
        The time at which the order is picked up from the restaurant
    delivery_time: int
        The time at which the order is delivered.
    courier_ID: str
        The string ID of the courier delivering the order.

    """

    def __init__(self, order: Order) -> None:
        """
        Create a new delivery instance.

        Parameters
        ----------
        order : Order
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.order_ID = order.code
        self.placement_time = order.placement_time
        self.ready_time = order.ready_time
        self.pickup_time = order.pickup_time
        self.delivery_time = order.delivery_time
        self.courier_ID = order.delivery_courier.code

    def __str__(self) -> str:
        """
        Return the text description of the delivery.

        The format of the description is as follows:
        Order_ID placement_time ready_time pickup_time delivery_time courier_ID

        Returns
        -------
        str
            The text description of the delivery.

        """
        return f'{self.order_ID} {self.placement_time} {self.ready_time} {self.pickup_time} {self.delivery_time} {self.courier_ID}'


class Dispatch:
    """A class holding information about courier dispatches.

    A courier dispatch represents a single movement of a courier from one
    location, one of home, restaurant or order dropoff location, to another
    location, either a restaurant or order dropoff location.

    Attributes
    ----------
    courier: Courier
        The courier being dispatched.
    departure_time: int
        The time at which the courier begins its dispatch.
    departure_location: Location
        The location at which the courier begins its dispatchment.
    arrival_location: Location
        The location at which the courier ends its dispatchment.

    """

    def __init__(self, courier: Courier, departure_time: int,
                 departure_location: Location, arrival_location: Location) -> None:
        """Create a new courier dispatch order.

        Parameters
        ----------
        courier : Courier
            The courier being dispatched.
        departure_time : int
            The time the courier begins its dispatchment.
        departure_location : Location
            The location the courier begins its dispatchment.
        arrival_location : Location
            The location the courier ends its dispatchment.

        Returns
        -------
        None.

        """
        self.courier = courier.code
        self.departure_time = departure_time
        self.departure_location = departure_location
        self.arrival_location = arrival_location

    def __str__(self) -> str:
        """
        Return a text description of the courier dispatchment.

        The format of the text description is:

            Courier_ID departure_time origin destination

        origin is '0' if the location of the origin is the courier's on
        location.

        Returns
        -------
        str
            The text description of the dispatch.

        """
        departure_location = 0 if type(self.departure_location) == Group \
            else self.departure_location.code
        return f'{self.courier} {self.departure_time} {departure_location} {self.arrival_location.code}'


assignment_information: [Assignment] = list()
order_delivery_information: [Delivery] = list()
courier_dispatch_information: [Dispatch] = list()


def display_solution() -> None:
    """
    Display the text version of the solution.

    Returns
    -------
    None.

    """
    print("assignment time, pickup time, Courier ID, Order ID, ..., Order ID")
    for assignment in assignment_information:
        print(assignment)
    print("\nOrder ID, placement time, ready time, pickup time, delivery \
time, courier ID")
    for delivery in order_delivery_information:
        print(delivery)
    print("\nCourier ID, departure time, origin, destination")
    for dispatch in courier_dispatch_information:
        print(dispatch)


class SolutionReader:
    """A class holding information to do with reading the model solution."""

    def parse_solution(fragments: [Fragment]) -> None:
        """
        Parse the given solution to make it readable.

        Parameters
        ----------
        fragments : [Fragment]
            A list of all activated fragments in the given model.

        Returns
        -------
        None.

        """
        for fragment in fragments:
            if len(fragment.order_list) > 0:
                assignment_information.append(Assignment(fragment))

            courier = fragment.group.couriers[0]
            if len(fragment.order_list) == 0:
                # In this case, the fragment only contains one courier
                # dispatch, and no order deliveries. If the departure location
                # and the arrival location are the same, then it is a waiting
                # fragment and we will ignore it. Otherwise, add the dispatch.
                if fragment.arrival_location == fragment.departure_location:
                    continue

                departure_time = fragment.departure_time
                departure_location = fragment.departure_location
                arrival_location = fragment.arrival_location
                courier_dispatch_information.append(Dispatch(
                    courier, departure_time, departure_location,
                    arrival_location))

            else:
                # The fragment contains at least one order, and so we need
                # to go through every order pair in the sequence, adding it
                # to the list of courier dispatches, and updating the new
                # information on the order class.
                departure_time = fragment.departure_time

                # Add the dispatch of the courier heading to the first order.
                courier_dispatch_information.append(
                    Dispatch(courier, departure_time,
                             fragment.departure_location,
                             fragment.order_list[0]))

                # Add a dispatch for every order. For the last order,
                # arrival_location is the next restaurant, unless it is the
                # courier's depot, in which case we ignore the final dispatch
                departure_time += fragment.departure_location.get_time_to(fragment.order_list[0]) + (Data.PICKUP_SERVICE_TIME + Data.DROPOFF_SERVICE_TIME) / 2

                # Iterate through the dispatches in the fragment
                for i in range(len(fragment.order_list) - 1):
                    # Store the delivery information for each order
                    order1 = fragment.order_list[i]
                    order1.delivery_time = departure_time
                    order1.pickup_time = fragment.departure_time
                    order1.delivery_courier = courier
                    order_delivery_information.append(Delivery(order1))

                    # Add a dispatch between the next two orders
                    order2 = fragment.order_list[i+1]
                    courier_dispatch_information.append(Dispatch(
                        courier, departure_time, order1, order2))

                    # Update the current time
                    departure_time += order1.get_time_to(order2) + Data.DROPOFF_SERVICE_TIME

                # Add a final dispatch to go from the last order to the
                # fragment's arrival_location, and add a delivery for the order
                last_order = fragment.order_list[-1]
                last_order.delivery_time = departure_time
                last_order.pickup_time = fragment.departure_time
                last_order.delivery_courier = courier
                order_delivery_information.append(Delivery(last_order))

                # Only add a dispatch for the final location if going to
                # another restaurant
                if type(fragment.arrival_location) == Restaurant:
                    courier_dispatch_information.append(Dispatch(
                        courier, departure_time, last_order,
                        fragment.arrival_location))

        # Sort the lists in the appropriate order
        # assignment_information should be sorted first by courier, then by
        # pickup_time
        # Sort assignments by pickup time
        assignment_information.sort(
            key=lambda assignment: assignment.pickup_time)
        # Sort assignments by courier ID
        assignment_information.sort(
            key=lambda assignment: assignment.courier_ID)
        # Sort deliveries by order ID
        order_delivery_information.sort(
            key=lambda delivery: delivery.order_ID)
        # Sort dispatches by courier ID
        courier_dispatch_information.sort(
            key=lambda dispatch: dispatch.courier)


class SolutionWriter:
    """A class to do with writing the solution to files."""

    def write_file(name: str, header: str, information: list()) -> None:
        """
        Create a new file and write in the given information.

        The function takes a given file name and header, and writes into the
        file given information. It does not write this information in the
        directory it is in, but rather in a folder called 'solution_files'.

        Parameters
        ----------
        name : str
            The name of the file to write.
        header : str
            The header line to write into the file.
        information : list()
            The list of information to write into the file.

        Returns
        -------
        None.

        """
        # Add the folder to the file name
        name = f'solution_files/{name}'
        with open(name, mode='w') as f:
            f.write(f'{header}\n')
            for line in information:
                # write each line in the information separately
                f.write(f'{line}\n')
            f.close()

    def write_solution() -> None:
        """
        Write three files, encompassing the entire model solution.

        Returns
        -------
        None.

        """
        # Write the assignment information
        SolutionWriter.write_file('assignment_solution_info.txt', 'assignment time, pickup time, Courier ID, Order ID, ..., Order ID', assignment_information)
        # Write the delivery information
        SolutionWriter.write_file('orders_solution_info.txt', 'Order ID, placement time, ready time, pickup time, delivery time, courier ID', order_delivery_information)
        # Write the dispatch information
        SolutionWriter.write_file('couriers_solution_info.txt', 'Courier ID, departure time, origin, destination', courier_dispatch_information)

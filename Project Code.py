"""
The main script running the program.

Created on Fri Feb 18 14:12:40 2022

@author: Tristan
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

The script follows a rough outline:
1. Import data
2. Compute sequences
3. Compute arcs
4. Create nodes
5. Create fragments
6. Create variables
7. Add objective function
8. Add constraints
9. Add callback logic
10. Optimise model
11. Convert solution to a readable format
12. Verify that the solution is in fact, valid

The following is a quick list of existing and planned functionality within the script.
# Done: Import data
# Done: Compute possible sequences
# Done: Create untimed arcs (courier/sequence/restaurant triples)
# Done: Create nodes per courier
# Done: Create timed arcs
# Done: Add model variables
# Done: Add objective function
# Done: Add model constraints
# Done: Solve model
# Done: Add cost penalty for not delivering orders
# Done: Add grouping couriers logic
# Done: Add increased discretisation logic

# Done: Have a node for each order ready time.
# Done: Create a map: If arrive at a restaruant at time t, snap to time t' (if active order, round down, otherwise, round up)
# Done: Strengthen constraint 2. Add new variable for if courier starts (Z_c) for example, then possibly sum constraint 2 with constraint 3.
# Done: Add a variable for courier starting, and change constraint 7 to take into account individual couriers, even in the group stage.

# Done: Add valid inequality logic
# Done: Remove variables with larger reduced cost.
# Done: Add only violated valid inequalities in loop
# Done: Add estimated optimal gurobi parameters
# TODO: In callback, if only add optimality cuts, check if solution used by gurobi, and if not, suggest it
# TODO: In callback, if add feasibility cuts, make a trivial solution by removing illegal orders
# TODO: In callback, if add feasibility cuts, solve mini IP with added fragments for group
# TODO: Find out why removing the bounds on any variable messes up objective value

Parameters
----------
consider_objective: Bool

    False for feasibility, True otherwise.

    This parameter affects the optimisation part of the model. If set to True, the model will solve
    the problem until it comes to an optimal solution to do with costs. If set to False, the model
    will only care about feasibility of the problem.

    NOTE: If set to False, parameter cost_penalty_active must be set to False. Otherwise the model
    will return a null result, that is, the model will not deliver any of the orders.

cost_penalty_active: Bool

    True to allow quick estimation of valid solutions, False otherwise.

    This parameter tells the model if it can consider solutions to the model where not all orders
    are delivered, at a large cost per order, or if it can only consider solutions where all orders
    are delivered.

    NOTE: If set to True, parameter consider_objective must be set to True. Otherwise, the model
    will return a null result, that is, the model will not deliver any of the orders.

courier_range_end: int

    The numeric ID of the first courier to include in the model.

    If parameter reduce_couriers is set to True, then the model will only include couriers whose
    numeric IDs fall within a run of integers. This run is bounded by parameters
    courier_range_start and courier_range_end to denote the start and end of this run, respectively.

    For example, if courier_range_start = 1 and courier_range_end = 5, then the model will consider
    couriers c1, c2, c3, c4 and c5, but will ignore all other couriers when running.

    NOTE: If parameter reduce_couriers = False, this parameter will have no effect.

courier_range_start: int

    The numeric ID of the first courier to include in the model.

    If parameter reduce_couriers is set to True, then the model will only include couriers whose
    numeric IDs fall within a run of integers. This run is bounded by parameters
    courier_range_start and courier_range_end to denote the start and end of this run, respectively.

    For example, if courier_range_start = 1 and courier_range_end = 5, then the model will consider
    couriers c1, c2, c3, c4 and c5, but will ignore all other couriers when running.

    NOTE: If parameter reduce_couriers = False, this parameter will have no effect.

couriers_to_avoid: [int]

    A list of numeric IDs of couriers to avoid.

    If parameter reduce_couriers is set to True, then the model will not include any orders whose
    numeric IDs fall within couriers_to_avoid.

    For example, if couriers_to_avoid = [1, 3, 5, 7, 9], then the model will not include orders c1,
    c3, c5, c7 or c9.

    Note: This parameter takes precedence over parameters courier_range_start and courier_range_end.
    For example, if courier_range_start = 1, courier_range_end = 3 and couriers_to_avoid = [2, 3,
    4], then the model will only include courier c1.

    Note: If parameter reduce_couriers = False, this parameter will have no effect.

group_couriers: Bool

    True to group couriers by their off time, false otherwise.

    This parameter determines whether the couriers are considered in terms of groups linked by
    their off times for the purposes of arc, node and fragment generation, or if each courier is
    considered its own group.

grubhub_instance: str

    A string denoting the instance being modelled.

    An instance string consists of five separate parts. The string takes the format of seed + size
    + speed + structure + prep. For example, if seed = 0, size = o50, speed = t75, structure = s1
    and prep = p100, then grubhub_instance = 0o50t75s1p100.

    seed is an integer ranging from 0 through to 9.

    size is either o100, denoting it contains the entire order set, o50, denoting it contains half
    of the order set, or r50, denoting it contains half of the restaurant set.

    speed is either t100, representing regular relative distances, or t75, representing close
    relative distances between locations.

    structure is either s1, representing days where couriers have similar shift durations and start
    and end times, or s2, representing days where couriers have arbitrary shift durations and start
    and end times.

    prep is either 100, representing low urgency orders or orders with short preparation times, or
    125, representing high urgency orders or orders with long preparation times.

order_range_end: int

    The numeric ID of the last order to include in the model.

    If parameter reduce_orders is set to True, then the model will only include orders whose
    numeric IDs fall within a run of integers. This run is bounded by parameters order_range_start
    and order_range_end to denote the start and end of this run, respectively.

    For example, if order_range_start = 1 and order_range_end = 5, then the model will consider
    orders o1, o2, o3, o4 and o5, but will ignore all other orders when running.

    NOTE: If parameter reduce_orders = False, this parameter will have no effect.

order_range_start: int

    The numeric ID of the first order to include in the model.

    If parameter reduce_orders is set to True, then the model will only include orders whose
    numeric IDs fall within a run of integers. This run is bounded by parameters order_range_start
    and order_range_end to denote the start and end of this run, respectively.

    For example, if order_range_start = 1 and order_range_end = 5, then the model will consider
    orders o1, o2, o3, o4 and o5, but will ignore all other orders when running.

    NOTE: If parameter reduce_orders = False, this parameter will have no effect.

orders_to_avoid: [int]

    A list of numeric IDs of orders to avoid.

    If parameter reduce_orders is set to True, then the model will not include any orders whose
    numeric IDs fall within orders_to_avoid.

    For example, if orders_to_avoid = [1, 3, 5, 7, 9], then the model will not include orders o1,
    o3, o5, o7 or o9.

    Note: This parameter takes precedence over parameters order_range_start and order_range_end.
    For example, if order_range_start = 1, order_range_end = 3 and orders_to_avoid = [2, 3, 4],
    then the model will only include order o1.

    Note: If parameter reduce_orders = False, this parameter will have no effect.

reduce_couriers: Bool

    True to ignore some of the couriers in the model, False otherwise.

reduce_orders: Bool

    True to ignore some of the orders in the model, False otherwise.

time_discretisation: int

    The amount of time between consecutive nodes with the same restaurant-group pair.

    This parameter affects the amount of time that passes between two nodes with the same
    restaurant and courier group coordinates. The minimum value is 1, and the estimated best value
    for computing time is equal to (2 + Data.PICKUP_SERVICE_TIME + Data.DROPOFF_SERVICE_TIME). Any
    values higher than the suggested value is at risk of having courier loops introduced in the
    optimisation part that need to be removed in the callback part.

"""

import time
import math
import itertools

from gurobipy import Model, quicksum, GRB
from classes import (
    Data,
    Courier,
    Order,
    Restaurant,
    Group,
    Sequence,
    Arc,
    Node,
    Fragment,
)
from typing import Set, Dict, Tuple, FrozenSet
from network import SubNetwork
from untimed_fragments_mdrp import UntimedFragmentsMDRP, ArcModelNoCourierAssignments

# Epoch time from program start
program_start_time = time.time()
grubhub_instance = "0o50t75s1p100"
# The directory containing the instance data
file_directory = "MealDeliveryRoutingGithub/public_instances/" + grubhub_instance + "/"

node_at_order_times = True
# TODO: implement this switch
time_discretisation = 10

reduce_orders = True
order_range_start = 1
order_range_end = 50
orders_to_avoid = set()

reduce_couriers = True
courier_range_start = 1  # TODO: Implement this functionality
courier_range_end = 61
couriers_to_avoid = list(
    (1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
     12, 13, 14, 16, 17, 19, 20,
     21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
     31, 32, 34, 35, 36, 37, 38, 40,
     41, 42, 43, 44, 46, 47, 48, 50,
     52, 53, 54, 55, 56, 58, 59, 60,
     61)
)
# Necessary couriers: 11, 15, 18, 33, 39, 45, 49, 51, 57

consider_objective = True
cost_penalty_active = True

group_by_off_time = True

""" ========== Valid Inequalities ========== """
add_valid_inequality_to_model = False
add_valid_inequality_after_LP = False
add_valid_inequality_to_callback = True

suggest_and_repair_solutions = True

# Output
summary_output = True
log_find_and_suggest_solutions = True
log_constraint_additions = True




def get_program_run_time() -> int:
    """Return the duration of the program as an integer."""
    return math.ceil(time.time() - program_start_time)


with open(file_directory + "instance_parameters.txt") as instance_parameters:
    """Read and save the instance parameters."""
    instance_parameters.readline().strip()
    parameters = instance_parameters.readline().strip().split("\t")
    Data.TRAVEL_SPEED = int(parameters[0])  # metres per minute
    Data.PICKUP_SERVICE_TIME = int(parameters[1])  # minutes
    Data.DROPOFF_SERVICE_TIME = int(parameters[2])  # minutes
    Data.TARGET_CLICK_TO_DOOR = int(parameters[3])  # minutes
    Data.MAX_CLICK_TO_DOOR = int(parameters[4])  # minutes
    Data.PAY_PER_DELIVERY = int(parameters[5])  # dollars
    Data.MIN_PAY_PER_HOUR = int(parameters[6])  # dollars

with open(file_directory + "couriers.txt") as file:
    lines = file.read().splitlines()[1:]
    for line in lines:
        data = line.split("\t")
        courier_number = int(data[0][1:])
        if not reduce_couriers or courier_number not in couriers_to_avoid:
            Courier(data[0], int(data[1]), int(data[2]), int(data[3]), int(data[4]))
print(f"{len(Courier.couriers)} total couriers")

with open(file_directory + "restaurants.txt") as file:
    lines = file.read().splitlines()[1:]
    for line in lines:
        data = line.split("\t")
        Restaurant(data[0], int(data[1]), int(data[2]))
print(f"{len(Restaurant.restaurants)} total restaurants")

with open(file_directory + "orders.txt") as file:
    lines = file.read().splitlines()[1:]
    for line in lines:
        data = line.split("\t")
        if reduce_orders:
            order_number = int(data[0][1:])
            if order_number >= order_range_start and order_number <= order_range_end:
                if order_number not in orders_to_avoid:
                    Order(
                        data[0],
                        int(data[1]),
                        int(data[2]),
                        int(data[3]),
                        data[4],
                        int(data[5]),
                    )
        else:
            Order(
                data[0], int(data[1]), int(data[2]), int(data[3]), data[4], int(data[5])
            )
print(f"{len(Order.orders)} total orders")
print(f"Data import complete at t = {get_program_run_time()}.\n")

print("Creating courier groups.")
Group.group_couriers(group_by_off_time)
print(f"Created {len(Group.groups)} groups at t = {get_program_run_time()}.\n")

print("Creating possible sequences.")
# Create sequences delivering orders
for restaurant in Restaurant.restaurants:
    if len(restaurant.orders) > 0:
        work_area = [Sequence([], restaurant)]
        while len(work_area) > 0:
            sequence = work_area[0]
            work_area.remove(sequence)
            for order in restaurant.orders:
                if order not in sequence.order_list:
                    new_sequence = sequence.add_order(order)
                    if (
                            new_sequence.latest_departure_time
                            >= new_sequence.earliest_departure_time
                    ):
                        work_area.append(new_sequence)
print("Sequence creation complete.")
print(
    f"Created {str(len(Sequence.sequences))} sequences at t = {get_program_run_time()}.\n"
)

print("Creating possible arcs.")
# Iterate through all non-empty sequences to create arcs
# For each non-empty sequence, iterate through all groups that can arrive in time to service the sequence
# Create an exit arc for that (group, sequence) pair
# Iterate through all restaurants
# If a restaurant contains an order that can be delivered after sequence, create an arc
for sequence in Sequence.sequences:
    if sequence.order_list != []:  # Create delivery sequences
        for group in Group.groups:
            earliest_departure_time = group.get_earliest_arrival_at(
                sequence.departure_location
            )
            if (
                    group.off_time >= sequence.earliest_departure_time
                    and earliest_departure_time <= sequence.latest_departure_time
            ):
                Arc(group, sequence, group)  # Creation of an exit arc
                for restaurant in Restaurant.restaurants:
                    earliest_arrival_time = (
                            earliest_departure_time
                            + sequence.travel_time
                            + restaurant.get_time_to(sequence.order_list[-1])
                            + (Data.DROPOFF_SERVICE_TIME + Data.PICKUP_SERVICE_TIME) / 2
                    )
                    deliverable_orders = set(
                        order
                        for order in restaurant.orders
                        if order not in sequence.order_list
                        and order.earliest_departure_time <= group.off_time
                        and order.latest_departure_time >= earliest_arrival_time
                    )
                    if len(deliverable_orders) > 0:
                        Arc(group, sequence, restaurant)
    else:  # Create wait arcs
        assert type(sequence.departure_location) == Restaurant
        arrival_restaurant: Restaurant
        arrival_restaurant = sequence.departure_location
        for group in Group.groups:
            earliest_arrival_time = group.get_earliest_arrival_at(arrival_restaurant)
            deliverable_orders = set(
                order
                for order in arrival_restaurant.orders
                if order.latest_departure_time >= earliest_arrival_time
                if order.earliest_departure_time <= group.off_time
            )
            if len(deliverable_orders) > 0:
                Arc(group, sequence, arrival_restaurant)
print("Created main, exit and wait arcs.")
for courier in Courier.couriers:  # entry arcs
    for restaurant in Restaurant.restaurants:
        earliest_arrival_time = courier.get_earliest_arrival_at(restaurant)
        deliverable_orders = set(
            order
            for order in restaurant.orders
            if order.earliest_departure_time <= courier.off_time
            and order.latest_departure_time >= earliest_arrival_time
        )
        if len(deliverable_orders) > 0:
            Arc(courier.get_group(), Sequence([], courier), restaurant)
print("Created entry arcs.")
print(f"Created {str(len(Arc.arcs))} arcs at t = {get_program_run_time()}.\n")

print("Creating model nodes.")
# Create a node for each (group, restaurant) pair at every earliest_departure_time of orders at that restaurant.
# Create a node for each courier at their on time.
# Create a node for each group at the off time.
# Node start time = max(group arrival, earliest deliverable order ready time)
# Node end time = min(group off time, latest deliverable order latest departure time)
for group in Group.groups:
    off_time = group.off_time
    for restaurant in Restaurant.restaurants:
        earliest_arrival_time = group.get_earliest_arrival_at(restaurant)
        deliverable_orders = set(
            order
            for order in restaurant.orders
            if order.earliest_departure_time <= off_time
            and order.latest_departure_time >= earliest_arrival_time
        )
        if len(deliverable_orders) > 0:
            node_start_time = int(
                max(
                    earliest_arrival_time,
                    min(order.earliest_departure_time for order in deliverable_orders),
                )
            )
            node_end_time = int(
                min(
                    off_time,
                    max(order.latest_departure_time for order in deliverable_orders),
                )
            )
            Node(group, restaurant, node_start_time)
            Node(group, restaurant, node_end_time)
            for order in deliverable_orders:
                if order.earliest_departure_time > node_start_time:
                    Node(group, restaurant, order.earliest_departure_time)
print("Restaurant nodes complete.")
for group in Group.groups:
    Node(group, group, group.off_time)
    for courier in group.couriers:
        Node(group, courier, courier.on_time)
print("Depot nodes complete.")
print(f"Created {str(len(Node.nodes))} nodes at t = {get_program_run_time()}.\n")

print("Creating path fragments.")
# Creation of path fragments that travel between locations. That is, deliver
# orders, or start at a depot and end at a restaurant
for arc in Arc.arcs:
    departure_nodes = Node.nodes_by_group_and_location[
        (arc.group, arc.departure_location)
    ]
    for node in departure_nodes:
        if (
                node.time >= arc.earliest_departure_time
                and node.time <= arc.latest_departure_time
        ):
            Fragment(node, arc)
print("Created path fragments.")

# # Creation of path fragments from depot to restaurants
# # If a valid depot -> restaurant fragment exists, then there is at least one node located at (group, restaurant)
# # Iterate through all (group, restaurant) pairs for nodes. Add a fragment from the depot to that (group, restaurant)
# for (group, restaurant) in Node.nodes_by_group_and_location:
#     node = Node.node_by_components[(group, group, 0)]
#     arc = Arc(group, Sequence([], 0, group, 0, 840, False), restaurant, False)
#     Fragment(node, arc, False)
# print('Created entry path fragments.')

# # Creation of 'waiting' path fragments.
# # For each (group, location) pair in the model, for every time associated
# # with that pair, create a path fragment from that node to the next one
# # chronologically.
# for (group, location) in Node.nodes_by_group_and_location:
#     arc = Arc(group, Sequence([], 0, location, time_discretisation, 840), location)
#     for node in Node.nodes_by_group_and_location[(group, location)]:
#         Fragment(node, arc, False)
# print('Created waiting path fragments.')
Fragment.sort()
print(
    f"Created {str(len(Fragment.fragments))} fragments at t = {get_program_run_time()}.\n"
)

print("Constructing model.")
mdrp = Model("Meal Delivery Routing Problem")
mdrp.setParam("LazyConstraints", 1)
mdrp.setParam("Method", 2)
mdrp.setParam("MIPGAP", 0.1)
# Tuning suggested parameters:
# mdrp.setParam('DegenMoves', 1)
# mdrp.setParam('Heuristics', 0.001)
# mdrp.setParam('MIPFocus', 1)
# mdrp.setParam('PrePasses', 1)

print("Creating variables.")
fragment_variables = {fragment: mdrp.addVar() for fragment in Fragment.fragments}
order_variables = {order: mdrp.addVar() for order in Order.orders}
payment_variables = {group: mdrp.addVar() for group in Group.groups}
courier_variables = {courier: mdrp.addVar(ub=1) for courier in Courier.couriers}

if consider_objective:
    print("Defining objective.")
    if cost_penalty_active:
        mdrp.setObjective(
            quicksum(payment_variables[group] for group in Group.groups)
            + 10000 * quicksum(1 - order_variables[order] for order in Order.orders)
        )
    else:
        mdrp.setObjective(quicksum(payment_variables[group] for group in Group.groups))

print("Creating constraints.")
if consider_objective:
    pay_for_deliveries = {
        group: mdrp.addConstr(
            payment_variables[group]
            >= quicksum(
                fragment_variables[fragment] * len(fragment.order_list) * Data.PAY_PER_DELIVERY
                for fragment in Fragment.fragments_by_group[group]
            )
            + quicksum(
                (1 - courier_variables[courier])
                * Data.MIN_PAY_PER_HOUR
                * (courier.off_time - courier.on_time)
                / 60
                for courier in group.couriers
            )
        )
        for group in Group.groups
    }

    pay_for_time = {
        group: mdrp.addConstr(
            payment_variables[group] >= Data.MIN_PAY_PER_HOUR * group.get_total_on_time() / 60
        )
        for group in Group.groups
    }

deliver_all_orders = {
    order: mdrp.addConstr(
        quicksum(fragment_variables[fragment] for fragment in Fragment.fragments_by_order[order])
        == order_variables[order]
    )
    for order in Order.orders
}

if cost_penalty_active is False:
    all_orders_on = {
        order: mdrp.addConstr(order_variables[order] == 1) for order in Order.orders
    }
else:
    all_orders_on = {
        order: mdrp.addConstr(order_variables[order] <= 1) for order in Order.orders
    }

flow_in_equals_flow_out = {
    node: mdrp.addConstr(
        quicksum(
            fragment_variables[fragment] for fragment in Fragment.fragments_by_arrival_node[node]
        )
        == quicksum(
            fragment_variables[fragment]
            for fragment in Fragment.fragments_by_departure_node[node]
        )
    )
    for node in Node.nodes
    if type(node.location) == Restaurant
}

couriers_start_once = {
    courier: mdrp.addConstr(
        quicksum(
            fragment_variables[fragment]
            for fragment in Fragment.departure_fragments_by_courier[courier]
        )
        == courier_variables[courier]
    )
    for courier in Courier.couriers
}

if add_valid_inequality_to_model:
    predecessor_valid_inequalities = {arc: mdrp.addConstr(quicksum(
        fragment_variables[fragment]
        for predecessor in arc.get_pred()
        for fragment in Fragment.fragments_by_arc[predecessor]
    )
                                                          >= quicksum(
        fragment_variables[fragment] for fragment in Fragment.fragments_by_arc[arc]
    ), name=f'predecessor valid inequality for {arc}')
                                      for arc in Arc.arcs
                                      if type(arc.departure_location) != Courier
                                      }
    successor_valid_inequalities = {arc: mdrp.addConstr(quicksum(
        fragment_variables[fragment]
        for successor in arc.get_succ()
        for fragment in Fragment.fragments_by_arc[successor]
    )
                                                        >= quicksum(
        fragment_variables[fragment] for fragment in Fragment.fragments_by_arc[arc]
    ), name=f'successor valid inequality for {arc}')
                                    for arc in Arc.arcs
                                    if type(arc.arrival_location) != Group
                                    }
print(f"Model finished constructing at t = {get_program_run_time()}.\n")

"""
Add Valid Inequalities to the model after solving the relaxed LP.

The program will repeatedly solve the LP, then check the arcs corresponding
with all activated fragments to check that they have the correct number of
predecessors and successors activated. If too few predecessors or successors
are activated, then a valid inequality will be added to the model, and the LP
solved again.
"""

VI = set()
if add_valid_inequality_after_LP:
    mdrp.Params.outputflag = 0
    while True:
        mdrp.optimize()

        VI_added: int = 0

        activated_arcs = set()
        # Find all activated arcs
        for fragment in fragment_variables:
            if fragment_variables[fragment].x > 0.01:
                activated_arcs.add(fragment.arc)

        # Add Valid Inequalities for the activated arcs
        for arc in activated_arcs:
            # Find value for arc
            arc_fragments = set(fragment for fragment in Fragment.fragments_by_arc[arc])
            arc_value = sum(fragment_variables[fragment].x for fragment in arc_fragments)
            # Add Valid Inequalities (if broken) on predecessors
            if type(arc.departure_location) is not Courier:
                pred_fragments = set(
                    fragment
                    for pred in arc.get_pred()
                    for fragment in Fragment.fragments_by_arc[pred]
                )
                pred_values = sum(fragment_variables[fragment].x for fragment in pred_fragments)
                if pred_values < arc_value - 0.1:
                    predecessor_VI = mdrp.addConstr(
                        quicksum(fragment_variables[fragment] for fragment in pred_fragments)
                        >= quicksum(fragment_variables[fragment] for fragment in arc_fragments)
                    )
                    VI_added += 1
                    VI.add(predecessor_VI)
            # Add Valid Inequalities (if broken) on successors
            if type(arc.arrival_location) is not Group:
                succ_fragments = set(
                    fragment
                    for succ in arc.get_succ()
                    for fragment in Fragment.fragments_by_arc[succ]
                )
                succ_values = sum(fragment_variables[fragment].x for fragment in succ_fragments)
                if succ_values < arc_value - 0.1:
                    successor_VI = mdrp.addConstr(
                        quicksum(fragment_variables[fragment] for fragment in succ_fragments)
                        >= quicksum(fragment_variables[fragment] for fragment in arc_fragments)
                    )
                    VI_added += 1
                    VI.add(successor_VI)
        print(
            f"Added {VI_added} violated valid inequalities at t = {get_program_run_time()}, objective = {mdrp.ObjVal}.")
        if VI_added == 0:
            break
    print(f'Added {len(VI)} Valid Inequalities in total at t = {get_program_run_time()}.\n')
    mdrp.Params.outputflag = 1

# lpgap = 0.01*mdrp.objVal
# rcList = [v for v in fragments.values() if v.rc > lpgap]
# print(len(rcList), "set to 0")
# mdrp.addConstr(quicksum(rcList) == 0)

for fragment in fragment_variables:
    fragment_variables[fragment].vtype = GRB.BINARY
for order in order_variables:
    order_variables[order].vtype = GRB.BINARY
for courier in courier_variables:
    courier_variables[courier].vtype = GRB.BINARY

mdrp._best_solution_value = GRB.INFINITY
mdrp._solved_subproblems = dict()


def get_length_of_order_overlap(fragment: Fragment, orders: FrozenSet[Order]) -> int:
    """Return the number of orders in the order list that overlap with the given order set."""
    return len(set(fragment.order_list).intersection(orders))


def add_lazy_optimality_cut_on_orders(model: Model, group: Group, orders: FrozenSet[Order], objective: float) -> None:
    """Add a lazy optimality cut on the given orders, for the given model."""
    arcs = Arc.get_arcs_for_orders(list(orders), group)
    cut_fragments = Fragment.get_fragments_from_arcs(arcs)
    model.cbLazy(
        payment_variables[group] >= objective * (quicksum(fragment_variables[fragment]
                                                          * get_length_of_order_overlap(fragment, orders)
                                                          for fragment in cut_fragments)
                                                 - len(orders) + 1))


def add_lazy_optimality_cut_on_arcs(model: Model, group: Group, arcs: FrozenSet[Arc], objective: float) -> None:
    """Add a lazy optimality cut on the given arcs, for the given model."""
    cut_fragments = Fragment.get_fragments_from_arcs(list(arcs))
    model.cbLazy(
        payment_variables[group] >= objective * (
                    quicksum(fragment_variables[fragment] for fragment in cut_fragments) - len(arcs) + 1))


def add_lazy_optimality_cut_on_fragments(model: Model, group: Group, fragments: FrozenSet[Fragment], objective: float) -> None:
    """Add a lazy optimality cut on the given fragments, for the given model."""
    model.cbLazy(
        payment_variables[group] >= objective * (
                    quicksum(fragment_variables[fragment] for fragment in fragments) - len(fragments) + 1))


def add_lazy_feasibility_predecessor_cut_on_arcs(model: Model, infeasible_arcs: Set[Arc], activated_arcs: Set[Arc]) -> None:
    """
    Add a lazy feasibility cut on the given infeasible_arcs, for the given model.

    The sum of the infeasible infeasible_arcs must be less than or equal to the number of
    infeasible infeasible_arcs minus 1, plus the sum of currently non-activated, alternative
    predecessors.
    """
    fragments_for_cut = set(Fragment.get_fragments_from_arcs(list(infeasible_arcs)))
    arc_predecessors = Arc.get_pred_to_arcs(infeasible_arcs)
    alternative_predecessors = set(arc_predecessors) - activated_arcs
    alternative_fragments = Fragment.get_fragments_from_arcs(list(alternative_predecessors))
    model.cbLazy(quicksum(fragment_variables[fragment] for fragment in fragments_for_cut)
                 <= len(fragments_for_cut) - 1
                 + quicksum(fragment_variables[fragment] for fragment in alternative_fragments))

def add_lazy_feasibility_cut_on_orders(model: Model, group: Group, feasible_orders: Set[Order],
                                       infeasible_orders: Set[Order]) -> None:
    """Add a lazy feasibility cut on the given orders, for the given model."""
    for order in infeasible_orders:
        orders_for_cut = list(feasible_orders)
        orders_for_cut.append(order)
        arcs_for_cut = Arc.get_arcs_with_orders(orders_for_cut, group)
        fragments_for_cut = Fragment.get_fragments_from_arcs(arcs_for_cut)
        model.cbLazy(quicksum(fragment_variables[fragment]
                              * get_length_of_order_overlap(fragment, frozenset(orders_for_cut))
                              for fragment in fragments_for_cut)
                     <= len(orders_for_cut) - 1)


def callback(model: Model, where: int) -> None:
    """A callback to handle feasibility and optimality cuts for the TF Model"""
    """
    for each group:
        find all assigned orders
        if (group, orders) pair already considered:

            if group objective value is too low:
                add optimality cut on group for orders
                save retrieved solution for suggestion
            else:
                save existing solution for suggestion

        else:
            find all used arcs for group
            solve UFModel on arcs
            if necessary, add feasibility cut on arcs
            if necessary, add optimality cut on arcs

            solve UFModel on orders
            if reached optimality:
                if necessary, add feasibility cut(s) on orders
                if necessary, add optimality cut on orders
                save (group, orders) pair as considered
            elif reached time limit:
                if necessary, add feasibility cut(s) on orders
                if necessary, add optimality cut on arcs

            if improved or arcs infeasible, save calculated solution for suggestion

    tally up group totals
    if group totals + extra sufficiently small:
        suggest all fragments
        suggest all group totals
    """
    if where == GRB.Callback.MIPSOL:
        summary_string = f'{model.cbGet(GRB.Callback.MIPSOL_OBJ)} -> '
        if not summary_output:
            print('-' * 50)
            print(f'Checking found solution of value {model.cbGet(GRB.Callback.MIPSOL_OBJ)}')
        # saved_solution[group] = (objective_value, extra_costs, solution_fragments)
        saved_solution: Dict[Group: Tuple[float, float, Set[Fragment]]] = dict()

        for group in Group.groups:
            if not summary_output:
                print(f'\nChecking solution for {group} of value {model.cbGetSolution(payment_variables[group])}')
            # find all assigned orders
            group_fragments = Fragment.fragments_by_group[group]
            group_fragments_solutions = model.cbGetSolution(
                list(fragment_variables[fragment] for fragment in group_fragments))
            gurobi_solution_group_fragments: Set[Fragment] = set()
            gurobi_solution_group_orders: Set[Order] = set()
            gurobi_solution_group_arcs: Set[Arc] = set()
            gurobi_solution_objective = model.cbGetSolution(payment_variables[group])
            for i in range(len(group_fragments_solutions)):
                if group_fragments_solutions[i] > 0.9:
                    gurobi_solution_group_fragments.add(group_fragments[i])
                    gurobi_solution_group_arcs.add(group_fragments[i].arc)
                    for order in group_fragments[i].order_list:
                        gurobi_solution_group_orders.add(order)
            if not summary_output:
                print(f'Gurobi assigned orders: {gurobi_solution_group_orders}')

            # check if (group, orders) pair already considered
            if (group, frozenset(gurobi_solution_group_orders)) in mdrp._solved_subproblems:
                if not summary_output:
                    print(f'Already solved on {(group, gurobi_solution_group_orders)}, comparing with existing solution')
                # if group objective value too low, add cut and suggest saved solution
                # if group objective is fine, save existing solution
                subnetwork = model._solved_subproblems[(group, frozenset(gurobi_solution_group_orders))]
                best_possible_solution_value = subnetwork.get_best_objective()
                extra_costs = subnetwork.get_undelivered_orders_cost()
                if best_possible_solution_value + extra_costs > gurobi_solution_objective + 0.01:
                    if not summary_output:
                        print(f'Found solution too low, adding optimality cut of {best_possible_solution_value}')
                    add_lazy_optimality_cut_on_orders(model, group, frozenset(gurobi_solution_group_orders),
                                                      best_possible_solution_value)
                    saved_solution[group] = (best_possible_solution_value, extra_costs, subnetwork.get_best_fragments())
                else:
                    if not summary_output:
                        print(f'Gurobi solution fine')
                    saved_solution[group] = (gurobi_solution_objective, 0, gurobi_solution_group_fragments)

            else:
                # Solve UFModel on arcs, and add optimality/feasibility cut if necessary
                subproblem = UntimedFragmentsMDRP(group, list(gurobi_solution_group_arcs),
                                                  list(gurobi_solution_group_orders))
                subproblem.optimize()
                if subproblem.getAttr(GRB.Attr.Status) == GRB.OPTIMAL:
                    # Add an optimality cut if necessary, else save gurobi's solution
                    objective_value = subproblem.getAttr(GRB.Attr.ObjVal)
                    if objective_value > gurobi_solution_objective + 0.01:
                        # Add an optimality cut
                        if not summary_output:
                            print('Gurobi solution too low, adding optimality cut on arcs')
                        add_lazy_optimality_cut_on_arcs(model, group, frozenset(gurobi_solution_group_arcs),
                                                        objective_value)
                        saved_solution[group] = (objective_value, 0, subproblem.convert_to_timed_path_fragments())
                    else:
                        if not summary_output:
                            print('Gurobi solution optimal on arcs')
                        saved_solution[group] = (gurobi_solution_objective, 0, gurobi_solution_group_fragments)

                else:
                    # subproblem infeasible, add feasibility cut
                    assert subproblem.getAttr(GRB.Attr.Status) == GRB.INFEASIBLE
                    if not summary_output:
                        print('Gurobi solution infeasible, adding feasibility cut on arcs')
                    infeasible_arcs = subproblem.get_infeasible_arcs()
                    add_lazy_feasibility_predecessor_cut_on_arcs(model, set(infeasible_arcs), gurobi_solution_group_arcs)
                    saved_solution[group] = (0, 10000 * len(gurobi_solution_group_orders), set())

                # Solve UFModel on orders, and add optimality/feasibility cut if necessary
                if not summary_output:
                    print(f'Looking for better solutions on orders for {group}')
                submodel = UntimedFragmentsMDRP.get_arc_model(group, list(gurobi_solution_group_orders))
                if submodel is None:
                    group_arcs = Arc.get_arcs_for_orders(list(gurobi_solution_group_orders), group)
                    submodel = UntimedFragmentsMDRP(group, group_arcs, list(gurobi_solution_group_orders),
                                                    cost_penalty_active=True, time_limit=10, save_model=True)
                    if not summary_output:
                        print(f'Created new submodel on {gurobi_solution_group_orders}')
                else:
                    if not summary_output:
                        print(f'Retrieved existing submodel on {gurobi_solution_group_orders}')
                    submodel.setParam('TimeLimit', submodel.getParamInfo('TimeLimit')[2] + 10)
                submodel.optimize()

                if not summary_output:
                    if submodel.getAttr(GRB.Attr.Status) == GRB.OPTIMAL:
                        print('Found optimal solution on group')
                    elif submodel.getAttr(GRB.Attr.Status) == GRB.TIME_LIMIT:
                        print('Time limit reached')
                    else:
                        print(f'Status = {submodel.getAttr(GRB.Attr.Status)}')

                # Add feasibility cuts if necessary
                undelivered_orders = submodel.get_undelivered_orders()
                delivered_orders = set(gurobi_solution_group_orders)
                delivered_orders -= set(undelivered_orders)
                if len(undelivered_orders) > 0 and submodel.getAttr(GRB.Attr.Status) == GRB.OPTIMAL:
                    if not summary_output:
                        print(f'Not all orders deliverable, cannot deliver {undelivered_orders}')
                    # Add feasibility cut on all undelivered orders
                    add_lazy_feasibility_cut_on_orders(model, group, delivered_orders, set(undelivered_orders))

                # Suggest found solution
                extra_costs = 10000 * len(undelivered_orders)
                objective_value = submodel.getAttr(GRB.Attr.ObjVal) - extra_costs
                solution_fragments = submodel.convert_to_timed_path_fragments()
                if objective_value + extra_costs + 0.01 < saved_solution[group][0] + saved_solution[group][1]:
                    if not summary_output:
                        print(f'Found solution of {objective_value + extra_costs}, suggesting new solution')
                    saved_solution[group] = (objective_value, extra_costs, solution_fragments)

                    """
                    Can I get to the point where I need to add an optimality cut?
                    An optimality cut should only be added if the best solution is 
                    of greater value than the current solution. Currently we know 
                    that we have a valid solution, and the value for that. So by
                    definition, the best solution has a value of less than or equal
                    to the currently known solution. 
                    Having said that, if the Gurobi solution is infeasible, I just
                    don't deliver any orders, and any solution I find that delivers 
                    orders is better.
                    """

                    # Add cut on orders if necessary
                    if submodel.getAttr(GRB.Attr.Status) == GRB.OPTIMAL:
                        if not summary_output:
                            print(f'Added optimality cut on orders of value {objective_value}')
                        add_lazy_optimality_cut_on_orders(model, group, delivered_orders, objective_value)
                    else:
                        if not summary_output:
                            print(f'Added optimality cut on fragments of value {objective_value}')
                        add_lazy_optimality_cut_on_fragments(model, group, solution_fragments, objective_value)

                # Solved to optimality, save submodel
                if submodel.getAttr(GRB.Attr.Status) == GRB.OPTIMAL:
                    if not summary_output:
                        print(f'Solved submodel to optimality with solution of {objective_value}')
                    mdrp._solved_subproblems[(group, frozenset(gurobi_solution_group_orders))] \
                        = SubNetwork(group, frozenset(gurobi_solution_group_orders), objective_value, extra_costs,
                                     set(solution_fragments))

        # Check to see if saved solution better than current best solution
        total_solution_value = 0
        for group in Group.groups:
            total_solution_value += saved_solution[group][0] + saved_solution[group][1]
        for order in Order.orders:
            if model.cbGetSolution(order_variables[order]) < 0.1:
                total_solution_value += 10000
        summary_string += f'{total_solution_value}'
        if total_solution_value + 0.01 < mdrp._best_solution_value:
            summary_string += ' (saved)'
            if not summary_output:
                print(f'New solution of {total_solution_value} best found so far')
            mdrp._best_solution_value = total_solution_value
            mdrp._best_solution_fragments = set()
            for group in Group.groups:
                mdrp._best_solution_fragments = mdrp._best_solution_fragments.union(saved_solution[group][2])
            mdrp._best_solution_group_values = {group: saved_solution[group][0] for group in Group.groups}
        else:
            if not summary_output:
                print(f'New solution of {total_solution_value} not best solution')
        if not summary_output:
            print('-' * 50)
        if summary_output:
            print(get_program_run_time(), summary_string)

    # Suggest saved solution to gurobi
    if where == GRB.Callback.MIPNODE:
        if mdrp._best_solution_value + 0.01 < model.cbGet(GRB.Callback.MIPNODE_OBJBST):
            if not summary_output:
                print(f'\nSuggesting solution of {mdrp._best_solution_value} to Gurobi')
            activated_fragments = mdrp._best_solution_fragments
            deactivated_fragments = set(Fragment.fragments) - activated_fragments
            model.cbSetSolution(list(fragment_variables[fragment] for fragment in activated_fragments),
                                list(1 for _ in range(len(activated_fragments))))
            model.cbSetSolution(list(fragment_variables[fragment] for fragment in deactivated_fragments),
                                list(0 for _ in range(len(deactivated_fragments))))
            model.cbSetSolution(list(payment_variables[group] for group in mdrp._best_solution_group_values),
                                list(mdrp._best_solution_group_values[group] for group in
                                     mdrp._best_solution_group_values))
            objVal = model.cbUseSolution()
            if not summary_output:
                print(f'Gurobi found objective of {objVal}')
                if objVal > mdrp._best_solution_value + 0.01:
                    for fragment in activated_fragments:
                        print(fragment)
            else:
                print(f'Suggested {mdrp._best_solution_value}')


# mdrp.setParam('TuneTimeLimit', 36000)
# mdrp.tune()
mdrp.optimize(callback)
print(f"\nProgram finished running at t = {get_program_run_time()}")



def orders_per_courier():
    """
    Group orders into sets according to their delivery courier.

    Returns
    -------
    orders_per_courier : {Courier: {Order}}
        A dictionary containing the orders delivered by each courier.

    """
    orders_per_courier = {courier: set() for courier in Courier.couriers}
    for fragment in Fragment.fragments:
        if fragment_variables[fragment].x > 0.9:
            courier = fragment.courier
            for order in fragment.order_list:
                orders_per_courier[courier].add(order)
    return orders_per_courier

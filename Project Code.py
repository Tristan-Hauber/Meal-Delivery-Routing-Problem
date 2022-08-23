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
from untimed_fragments_mdrp import UntimedFragmentsMDRP

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
order_range_end = order_range_start + 74
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

add_valid_inequality_to_model = False
add_valid_inequality_after_LP = False
add_valid_inequality_to_callback = False

suggest_and_repair_solutions = False

# Logging options
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
# Tuning suggested parameters:
# mdrp.setParam('DegenMoves', 1)
# mdrp.setParam('Heuristics', 0.001)
# mdrp.setParam('MIPFocus', 1)
# mdrp.setParam('PrePasses', 1)

print("Creating variables.")
fragments = {fragment: mdrp.addVar(ub=1) for fragment in Fragment.fragments}
orders = {order: mdrp.addVar(ub=1) for order in Order.orders}
payments = {group: mdrp.addVar() for group in Group.groups}
couriers = {courier: mdrp.addVar(ub=1) for courier in Courier.couriers}

if consider_objective:
    print("Defining objective.")
    if cost_penalty_active:
        mdrp.setObjective(
            quicksum(payments[group] for group in Group.groups)
            + 10000 * quicksum(1 - orders[order] for order in Order.orders)
        )
    else:
        mdrp.setObjective(quicksum(payments[group] for group in Group.groups))

print("Creating constraints.")
if consider_objective:
    pay_for_deliveries = {
        group: mdrp.addConstr(
            payments[group]
            >= quicksum(
                fragments[fragment] * len(fragment.order_list) * Data.PAY_PER_DELIVERY
                for fragment in Fragment.fragments_by_group[group]
            )
            + quicksum(
                (1 - couriers[courier])
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
            payments[group] >= Data.MIN_PAY_PER_HOUR * group.get_total_on_time() / 60
        )
        for group in Group.groups
    }

deliver_all_orders = {
    order: mdrp.addConstr(
        quicksum(fragments[fragment] for fragment in Fragment.fragments_by_order[order])
        == orders[order]
    )
    for order in Order.orders
}

if cost_penalty_active is False:
    all_orders_on = {
        order: mdrp.addConstr(orders[order] == 1) for order in Order.orders
    }
else:
    all_orders_on = {
        order: mdrp.addConstr(orders[order] <= 1) for order in Order.orders
    }

flow_in_equals_flow_out = {
    node: mdrp.addConstr(
        quicksum(
            fragments[fragment] for fragment in Fragment.fragments_by_arrival_node[node]
        )
        == quicksum(
            fragments[fragment]
            for fragment in Fragment.fragments_by_departure_node[node]
        )
    )
    for node in Node.nodes
    if type(node.location) == Restaurant
}

couriers_start_once = {
    courier: mdrp.addConstr(
        quicksum(
            fragments[fragment]
            for fragment in Fragment.departure_fragments_by_courier[courier]
        )
        == couriers[courier]
    )
    for courier in Courier.couriers
}

if add_valid_inequality_to_model:
    valid_inequalities = {
        arc: mdrp.addConstr(
            quicksum(
                fragments[fragment]
                for pred in arc.get_pred()
                for fragment in Fragment.fragments_by_arc[pred]
            )
            >= quicksum(
                fragments[fragment] for fragment in Fragment.fragments_by_arc[arc]
            )
        )
        for arc in Arc.arcs
        if type(arc.departure_location) != Courier
    }
print(f"Model finished constructing at t = {get_program_run_time()}.\n")

while True:

    mdrp.optimize()

    VI_added: int = 0

    if add_valid_inequality_after_LP:
        activated_arcs = set()
        for fragment in fragments:
            if fragments[fragment].x > 0.01:
                activated_arcs.add(fragment.arc)
        for arc in activated_arcs:
            pred_fragments = set(
                fragment
                for pred in arc.get_pred()
                for fragment in Fragment.fragments_by_arc[pred]
            )
            succ_fragments = set(
                fragment
                for succ in arc.get_succ()
                for fragment in Fragment.fragments_by_arc[succ]
            )
            arc_fragments = set(fragment for fragment in Fragment.fragments_by_arc[arc])

            pred_values = sum(fragments[fragment].x for fragment in pred_fragments)
            succ_values = sum(fragments[fragment].x for fragment in succ_fragments)
            arc_value = sum(fragments[fragment].x for fragment in arc_fragments)

            if pred_values < arc_value:
                mdrp.addConstr(
                    quicksum(fragments[fragment] for fragment in pred_fragments)
                    >= quicksum(fragments[fragment] for fragment in arc_fragments)
                )
                VI_added += 1

            if succ_values > arc_value:
                mdrp.addConstr(
                    quicksum(fragments[fragment] for fragment in succ_fragments)
                    >= quicksum(fragments[fragment] for fragment in arc_fragments)
                )
                VI_added += 1
        print(f"Added {VI_added} violated valid inequalities.")

    if VI_added == 0:
        break

# lpgap = 0.01*mdrp.objVal
# rcList = [v for v in fragments.values() if v.rc > lpgap]
# print(len(rcList), "set to 0")
# mdrp.addConstr(quicksum(rcList) == 0)

for fragment in fragments:
    fragments[fragment].vtype = GRB.BINARY
for order in orders:
    orders[order].vtype = GRB.BINARY
for courier in couriers:
    couriers[courier].vtype = GRB.BINARY

mdrp._best_solution_value = GRB.INFINITY


def Callback(model, where):
    """
    Check provided solutions for legality.

    This callback provides one purpose, and that is to eliminate paths from the network. There are two non-exclusive occasions in which this might be needed; one is when the couriers are grouped by off time, and the other is when nodes are further apart than integers. In both cases, a sub-network is created from all activated arcs and solved, before adding cuts to the main model.

    If the couriers are grouped by off time, we are guaranteed to have a feasible solution, but we are not guaranteed that is optimal. In this case, we have a sub-network for each courier group, which is solved to optimality. Then an optimality cut is imposed onto the main problem to reflect the solution found in the sub-network.

    If the timing points on nodes are further apart than integers, we are not guaranteed to have a feasible solution. In this case, we need to create a sub-network per courier and solve it to determine if the solution is feasible. If it is not feasible, we add a feasibility cut to the main model.

    If the couriers are grouped by off time and the timings on nodes are further apart than integers, we are not guaranteed a feasible solution, and if we do get a feasible solution, the solution may not be optimal. In this case, we combine the methods from the individual cases, and solve a sub-network for each courier group for feasibility first and optimality second, adding feasibility and optimality cuts to the main model where needed.

    Parameters
    ----------
    model : Model
        The model the callback is being called on.
    where : int
        The current location in the solve.

    Returns
    -------
    None.

    """
    if where == GRB.Callback.MIPSOL:
        if log_find_and_suggest_solutions:
            print(
                f"Checking new incumbent solution with value {mdrp.cbGet(GRB.Callback.MIPSOL_OBJ)}"
            )
        suggest_solution = suggest_and_repair_solutions
        # Get all activated fragments
        group_arcs = {group: [] for group in Group.groups}
        group_orders = {group: [] for group in Group.groups}
        group_fragments = {group: [] for group in Group.groups}
        group_payments = {
            group: group.get_total_on_time() * Data.MIN_PAY_PER_HOUR / 60
            for group in Group.groups
        }

        for fragment in fragments:

            if mdrp.cbGetSolution(fragments[fragment]) < 0.1:
                continue

            # Add all activated fragments to the list, even waiting fragments
            group_fragments[fragment.group].append(fragment)

            # We don't care about any fragments that aren't moving between locations, that is, waiting fragments
            if (
                len(fragment.order_list) == 0
                and fragment.departure_location == fragment.arrival_location
            ):
                continue

            # Add activated arcs to the list
            group = fragment.group
            group_arcs[group].append(fragment.arc)
            for order in fragment.order_list:
                group_orders[group].append(order)

        # We have identified all activated fragments
        # Now we go through each courier group and build a sub-network
        for group in Group.groups:
            arcs = group_arcs[group]

            # If no activated arcs, then subnetwork is trivially feasible and optimal
            if len(arcs) == 0:
                continue

            """
            The following two for loops set up the predecessor and successor dictionaries for later use. These are useful as in order to solve the network, we need to keep track of predecessors and successors within the network.

            For arc pred to be a predecessor of arc succ, then pred must finish at the same location that succ begins, and the courier must be able to arrive at that location before they must leave for succ.

            For arc succ to be a successor of arc pred, then succ must begin at the same location that pred ends, and the courier must have time to finish pred before they have to leave for succ.

            From these definitions, we can see that an arc pred is a predecessor of an arc succ if and only if succ is a successor of pred.
            """
            # Dictionaries for storing predecessors and arcs.
            successors_of_arc = {arc: set() for arc in arcs if arc.arrival_location != arc.group}

            for arc1, arc2 in itertools.combinations(arcs, 2):
                # Check if arc2 is a successor of arc1
                if arc1.arrival_location != arc1.group and arc1.arrival_location == arc2.departure_location and arc1.earliest_departure_time + arc1.travel_time <= arc2.latest_departure_time:
                    # Add arc2 as a successor of arc1
                    successors_of_arc[arc1].add(arc2)
                # Check if arc1 is a successor of arc2
                if arc2.arrival_location != arc2.group and arc2.arrival_location == arc1.departure_location and arc2.earliest_departure_time + arc2.travel_time <= arc1.latest_departure_time:
                    # Add arc1 as a successor of arc2
                    successors_of_arc[arc2].add(arc1)

            """
            The rest of this function deals with the creation and solving of the sub-model. The mathematics of this sub-model are explained in detail in sections 3.1 and 3.2 of the write-up.

            We start off by constructing the variables. We have one variable for courier payments within the group, one variable keeping track of which arcs are used as successors or predecessors, one variable keeping track of when each arc is serviced, and one variable linking arcs with servicing couriers.

            The objective is to minimise the total payment to the couriers within the group.

            Then we have the constraints linking the variables together. Note: the equation numbers given in the code may not be correct.
            """
            ipe = Model("Illegal Path Elimination")
            ipe.setParam("OutputFlag", 0)
            courier_payments = {courier: ipe.addVar() for courier in group.couriers}
            successors = {(arc1, arc2): ipe.addVar(vtype=GRB.BINARY) for arc1 in successors_of_arc for arc2 in successors_of_arc[arc1]}
            timings = {arc: ipe.addVar() for arc in arcs}
            assignments = {
                (courier, arc): ipe.addVar(vtype=GRB.BINARY)
                for courier in group.couriers
                for arc in arcs
            }

            # Minimise cumulative courier payments (eq 8)
            ipe.setObjective(
                quicksum(courier_payments[courier] for courier in group.couriers)
            )

            # Each courier is paid for their time (eq 9)
            courier_time = {
                courier: ipe.addConstr(
                    courier_payments[courier]
                    >= (courier.off_time - courier.on_time) / 60 * Data.MIN_PAY_PER_HOUR
                )
                for courier in group.couriers
            }
            # Each courier is paid per delivery (eq 10)
            courier_deliveries = {
                courier: ipe.addConstr(
                    courier_payments[courier]
                    >= quicksum(
                        assignments[courier, arc]
                        * len(arc.order_list)
                        * Data.PAY_PER_DELIVERY
                        for arc in arcs
                    )
                )
                for courier in group.couriers
            }
            # All activated arcs are assigned to a courier (eq 11)
            arcs_serviced = {
                arc: ipe.addConstr(
                    quicksum(assignments[courier, arc] for courier in group.couriers)
                    == 1
                )
                for arc in arcs
            }
            # All non-exit arcs must have a successor (eq 12)
            have_succ = {arc1: ipe.addConstr(quicksum(successors[arc1, arc2] for arc2 in successors_of_arc[arc1]) == 1) for arc1 in arcs if arc1.arrival_location != arc1.group}
            # All successors must start late enough for the previous to finish (eq 15)
            succ_timings = {
                (arc1, arc2): ipe.addConstr(
                    timings[arc1] + arc1.travel_time
                    <= timings[arc2]
                    + (1 - successors[arc1, arc2])
                    * (
                        arc1.latest_departure_time
                        + arc1.travel_time
                        - arc2.earliest_departure_time
                    )
                )
                for (arc1, arc2) in successors
            }
            # If a courier services an arc, it must also service its successor (eq 16)
            cour_serv_succ = {
                (arc1, arc2, courier): ipe.addConstr(
                    assignments[courier, arc1] + successors[arc1, arc2] - 1
                    <= assignments[courier, arc2]
                )
                for courier in group.couriers
                for (arc1, arc2) in successors
            }
            # An arc must be serviced after it is ready (eq 13)
            begin_on_time = {
                arc: ipe.addConstr(timings[arc] >= arc.earliest_departure_time)
                for arc in arcs
            }
            # An arc must be serviced before it is too late (eq 14)
            end_on_time = {
                arc: ipe.addConstr(timings[arc] <= arc.latest_departure_time)
                for arc in arcs
            }

            # After constraint generation, solve submodel
            ipe.optimize()

            if ipe.getAttr(GRB.Attr.Status) == GRB.OPTIMAL:
                # Sub-model is feasible, move to optimality cuts
                group_payment = ipe.getObjective().getValue()
                group_payments[group] = group_payment
                # If optimal value greater than calculated value, add optimality cut
                if group_payment > mdrp.cbGetSolution(payments[group]):
                    # Equation 17
                    mdrp.cbLazy(
                        payments[group]
                        >= group_payment
                        * (
                            1
                            - len(arcs)
                            + quicksum(
                                fragments[fragment]
                                for arc in arcs
                                for fragment in Fragment.fragments_by_arc[arc]
                            )
                        )
                    )
                    if log_constraint_additions:
                        print(f"Added optimality constraint on {group}")

            else:
                # Sub-model is infeasible, move to feasibility cuts
                ipe.computeIIS()
                infeasible_arcs = set()
                # Not all groups are feasible

                # Add all infeasible arcs to our collection
                for arc in arcs_serviced:
                    if arcs_serviced[arc].IISConstr:
                        infeasible_arcs.add(arc)
                for arc in have_succ:
                    if have_succ[arc].IISConstr:
                        infeasible_arcs.add(arc)
                for (arc1, arc2) in succ_timings:
                    if succ_timings[(arc1, arc2)].IISConstr:
                        infeasible_arcs.add(arc1)
                        infeasible_arcs.add(arc2)
                for (arc1, arc2, cour) in cour_serv_succ:
                    if cour_serv_succ[(arc1, arc2, cour)].IISConstr:
                        infeasible_arcs.add(arc1)
                        infeasible_arcs.add(arc2)
                for arc in begin_on_time:
                    if begin_on_time[arc].IISConstr:
                        infeasible_arcs.add(arc)
                for arc in end_on_time:
                    if end_on_time[arc].IISConstr:
                        infeasible_arcs.add(arc)

                # Add a feasibility cut (eq 18)
                mdrp.cbLazy(
                    quicksum(
                        fragments[fragment]
                        for arc in infeasible_arcs
                        for fragment in Fragment.fragments_by_arc[arc]
                    )
                    <= len(infeasible_arcs)
                    - 1
                    + quicksum(
                        fragments[fragment]
                        for pred in Arc.get_pred_to_arcs(infeasible_arcs)
                        if pred not in arcs
                        for fragment in Fragment.fragments_by_arc[pred]
                    )
                )
                if log_constraint_additions:
                    print(f"Added feasibility cut on {group}")

                # Add a valid inequality cut
                if add_valid_inequality_to_callback:
                    VI_added: int = 0
                    for arc in infeasible_arcs:
                        mdrp.cbLazy(
                            quicksum(
                                fragments[fragment]
                                for pred in arc.get_pred()
                                for fragment in Fragment.fragments_by_arc[pred]
                            )
                            >= quicksum(
                                fragments[fragment]
                                for fragment in Fragment.fragments_by_arc[arc]
                            )
                        )
                        mdrp.cbLazy(
                            quicksum(
                                fragments[fragment]
                                for succ in arc.get_succ()
                                for fragment in Fragment.fragments_by_arc[succ]
                            )
                            >= quicksum(
                                fragments[fragment]
                                for fragment in Fragment.fragments_by_arc[arc]
                            )
                        )
                        VI_added += 2
                    if log_constraint_additions:
                        print(f"Added {VI_added} valid inequalities.")

                if not suggest_solution:
                    continue
                # Create a new solution for the group
                if log_find_and_suggest_solutions:
                    print(f'Searching for new solution for {group} with orders {group_orders[group]}.')
                submodel = UntimedFragmentsMDRP.get_ufmdrp(group, group_orders[group])
                # Check to see if already created submodel
                if submodel is None:
                    submodel = UntimedFragmentsMDRP(group, group.get_arcs(), group_orders[group], save_model=True)
                else:
                    if log_find_and_suggest_solutions:
                        print('Retrieving existing model.')
                    submodel.uf_mdrp.Params.time_limit += 5
                submodel.optimize()
                # Only proceed if we have at least one solution
                if submodel.uf_mdrp.getAttr('Status') == GRB.OPTIMAL:
                    if log_find_and_suggest_solutions:
                        print(f'Solution found for {group}.')
                    # Only add a constraint if the model solved close enough to optimality
                    # Add constraint on orders in the group
                    undelivered_orders = submodel.get_undelivered_orders()
                    if len(undelivered_orders) > 0:
                        delivered_orders = list(order for order in submodel.deliveries if order not in undelivered_orders)
                        abs_gap = submodel.uf_mdrp.ObjVal - submodel.uf_mdrp.ObjBound
                        worst_reduction = 10000 - Data.PAY_PER_DELIVERY
                        deliverable_orders = math.ceil(abs_gap / worst_reduction)
                        # TODO: Confirm that this maths checks out
                        if deliverable_orders < len(undelivered_orders):
                            for order_combination in itertools.combinations(undelivered_orders, deliverable_orders + 1):
                                invalid_order_set = set(delivered_orders)
                                invalid_order_set.union(order_combination)
                                invalid_fragment_set = set(fragment for fragment in Fragment.get_fragments_from_orders(invalid_order_set) if fragment.group == group)
                                mdrp.cbLazy(quicksum(fragments[fragment] * len(set(fragment.order_list).intersection(invalid_order_set)) for fragment in invalid_fragment_set) <= len(invalid_order_set) - 1)
                            if log_constraint_additions:
                                print(f'Limited {group} to {delivered_orders} and no group of {deliverable_orders + 1} orders from {undelivered_orders}.')
                    # Suggest a solution to Gurobi
                    if suggest_solution:
                        group_payments[group] = submodel.objVal
                        group_fragments[group] = submodel.convert_to_timed_path_fragments()
                        if log_find_and_suggest_solutions:
                            print(f'Found new solution for {group}.')
                else:
                    if log_find_and_suggest_solutions:
                        print(f'No solution found for {group}. Will not suggest a solution.')
                    suggest_solution = False

        """
        Suggest our constructed solution to Gurobi.

        If all groups are feasible, then we simply suggest the fragments that were activated at the beginning of the callback. If at least one group was infeasible, then we replace the fragments for that group with our new fragments.
        """
        # Calculate new solution value
        if suggest_solution:
            solution_value = 0
            for group in Group.groups:
                solution_value += group_payments[group]
            for order in orders:
                if mdrp.cbGetSolution(orders[order]) < 0.1:
                    solution_value += 10000
            # Compare to previous solution value
            if solution_value + 0.0001 < mdrp._best_solution_value:
                # Save activated fragments and delivered orders to suggest later
                mdrp._best_fragments = list()
                for group in group_fragments:
                    mdrp._best_fragments += group_fragments[group]
                mdrp._best_solution_value = solution_value
                if log_find_and_suggest_solutions:
                    print(
                        f"Saved solution of value {solution_value} to suggest later\n"
                    )
            elif log_find_and_suggest_solutions:
                print(f"Found solution of {solution_value} worse than current best solution\n")
        else:
            if log_find_and_suggest_solutions:
                print()

    """
    In the MIPNODE stage of solving the problem, we can suggest previous
    solutions that we have found. We check that our suggested solution is
    better than the currently best found solution, and then suggest it to
    Gurobi to implement into the path.
    """
    if (
        where == GRB.Callback.MIPNODE
        and mdrp._best_solution_value + 0.0001 < mdrp.cbGet(GRB.Callback.MIPNODE_OBJBST)
    ):
        # Give Gurobi values for all fragments according to our solution
        activated_fragment_variables = list()
        non_activated_fragment_variables = list()
        for fragment in fragments:
            if fragment in mdrp._best_fragments:
                activated_fragment_variables.append(fragments[fragment])
            else:
                non_activated_fragment_variables.append(fragments[fragment])
        mdrp.cbSetSolution(activated_fragment_variables, list(1 for _ in range(len(activated_fragment_variables))))
        mdrp.cbSetSolution(non_activated_fragment_variables, list(0 for _ in range(len(non_activated_fragment_variables))))
        objVal = mdrp.cbUseSolution()
        if log_find_and_suggest_solutions:
            print(f'Suggested value of {mdrp._best_solution_value}, Gurobi found {objVal}')
            print('Activated fragments are:')
            for fragment in mdrp._best_fragments:
                print(fragment)
            print()


# mdrp.setParam('TuneTimeLimit', 36000)
# mdrp.tune()
mdrp.optimize(Callback)
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
        if fragments[fragment].x > 0.9:
            courier = fragment.courier
            for order in fragment.order_list:
                orders_per_courier[courier].add(order)
    return orders_per_courier

"""
The main script running the program.

Created on Fri Feb 18 14:12:40 2022

@author: Tristan

The script follows a rough outline:
1. Import data
2. Create sets
3. Create variables
4. Add objective function
5. Add constraints
6. Define callback
7. Optimise model
8. Verify validity of solution

Parameters

grubhub_instance: str

reduce_orders: Bool
order_range_start: int
order_range_end: int
orders_to_avoid: int

reduce_couriers: Bool
courier_range_start: int
courier_range_end: int
couriers_to_avoid: int

maximum_sequence_length: int

node_at_order_times: Bool
time_discretisation: int

group_by_off_time: Bool

consider_objective: Bool
cost_penalty_active: Bool
cost_penalty_value: int

add_valid_inequality_to_model: Bool
add_valid_inequality_after_LP: Bool
add_valid_inequality_to_callback: Bool

submodel_version: int

suggest_solutions: Bool
repair_solutions: Bool
find_new_solutions: Bool

model_output_type = "Summary"
submodel_output_type = "None"
"""

import time

from gurobipy import Model, quicksum, GRB
from classes import *
from typing import Set, Dict, Tuple, FrozenSet, List
from network import SubNetwork
from untimed_fragments_mdrp import UntimedFragmentsMDRP
from uf_mdrp_courier_successors import UFMDRP2
from uf_mdrp_3 import UFMDRP3
import solution_checker

program_start_time = time.time()

""" ========== INSTANCE PARAMETERS ========== """
grubhub_instance = "0o50t75s1p100"
file_directory = "MealDeliveryRoutingGithub/public_instances/" + grubhub_instance + "/"
# Order Generation
reduce_orders = True
order_range_start = 1
order_range_end = 20
orders_to_avoid = set()
# Courier Generation
reduce_couriers = True
courier_range_start = 1
courier_range_end = 61
couriers_to_avoid = list(
    (1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
     12, 13, 14, 16, 17, 19, 20,
     21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
     31, 32, 34, 35, 36, 37, 38, 40,
     41, 42, 43, 44, 46, 47, 48, 50,
     52, 53, 54, 55, 56, 58, 59, 60,
     61)
)  # Necessary couriers: 11, 15, 18, 33, 39, 45, 49, 51, 57
# Sequence Generation
maximum_sequence_length = 2

""" ========== MODEL PARAMETERS ========== """
# Node Generation
node_at_order_times = True
time_discretisation = 1
# Group Generation
group_by_off_time = True
# Objective Generation
consider_objective = True
cost_penalty_active = True
cost_penalty_value = 10000
# Valid Inequalities
add_valid_inequality_to_model = False
add_valid_inequality_after_LP = False
add_valid_inequality_to_callback = False

""" ========== SUBMODEL PARAMETERS ========== """
# Submodel Formulation
submodel_version = 1
if submodel_version == 1:
    SubModel = UntimedFragmentsMDRP
elif submodel_version == 2:
    SubModel = UFMDRP2
elif submodel_version == 3:
    SubModel = UFMDRP3
else:
    assert False
# Solution Searching
suggest_solutions = True
repair_solutions = False
find_new_solutions = True
suggest_and_repair_solutions = True

""" ========== OUTPUT ========== """
model_output_type = "Summary"  # Full, Summary or None
submodel_output_type = "None"  # Full, Summary or None
summary_output = True
log_find_and_suggest_solutions = True
log_constraint_additions = True


def get_program_run_time() -> float:
    """Return the duration of the program as an integer."""
    return round(time.time() - program_start_time, 2)


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
        if reduce_couriers:
            courier_number = int(data[0][1:])
            if courier_number in couriers_to_avoid \
                    or courier_number < courier_range_start \
                    or courier_number > courier_range_end:
                continue
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
            if order_number in orders_to_avoid \
                    or order_number < order_range_start \
                    or order_number > order_range_end:
                continue
        Order(data[0], int(data[1]), int(data[2]), int(data[3]), data[4], int(data[5]))
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
                if order not in sequence.order_list and len(sequence.order_list) < maximum_sequence_length:
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
    if len(sequence.order_list) > 0:  # Create delivery sequences
        for courier_group in Group.groups:
            earliest_departure_time = courier_group.get_earliest_arrival_at(
                sequence.departure_location
            )
            if (
                    courier_group.off_time >= sequence.earliest_departure_time
                    and earliest_departure_time <= sequence.latest_departure_time
            ):
                Arc(courier_group, sequence, courier_group)  # Creation of an exit arc
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
                        and order.earliest_departure_time <= courier_group.off_time
                        and order.latest_departure_time >= earliest_arrival_time
                    )
                    if len(deliverable_orders) > 0:
                        Arc(courier_group, sequence, restaurant)
    else:  # Create wait arcs
        assert type(sequence.departure_location) == Restaurant
        arrival_restaurant: Restaurant
        arrival_restaurant = sequence.departure_location
        for courier_group in Group.groups:
            earliest_arrival_time = courier_group.get_earliest_arrival_at(arrival_restaurant)
            deliverable_orders = set(
                order
                for order in arrival_restaurant.orders
                if order.latest_departure_time >= earliest_arrival_time
                if order.earliest_departure_time <= courier_group.off_time
            )
            if len(deliverable_orders) > 0:
                Arc(courier_group, sequence, arrival_restaurant)
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
# Node start time = max(group arrival, earliest departure time of the earliest deliverable order)
# Node end time = min(group off time, latest departure time of the latest deliverable order)
for courier_group in Group.groups:
    off_time = courier_group.off_time
    for restaurant in Restaurant.restaurants:
        earliest_arrival_time = courier_group.get_earliest_arrival_at(restaurant)
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
            if node_at_order_times:
                Node(courier_group, restaurant, node_start_time)
                Node(courier_group, restaurant, node_end_time)
                for order in deliverable_orders:
                    if order.earliest_departure_time > node_start_time:
                        Node(courier_group, restaurant, order.earliest_departure_time)

            else:
                node_time = node_start_time
                while node_time <= node_end_time:
                    Node(courier_group, restaurant, node_time)
                    node_time += time_discretisation

print("Restaurant nodes complete.")
for courier_group in Group.groups:
    Node(courier_group, courier_group, courier_group.off_time)
    for courier in courier_group.couriers:
        Node(courier_group, courier, courier.on_time)
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
        if arc.earliest_departure_time <= node.time <= arc.latest_departure_time:
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
fragment_variables = {fragment: mdrp.addVar(lb=0, ub=GRB.INFINITY, obj=0, vtype=GRB.CONTINUOUS,
                                            name=f'{fragment} activation variable', column=None)
                      for fragment in Fragment.fragments}
order_variables = {order: mdrp.addVar(lb=0, ub=GRB.INFINITY, obj=0, vtype=GRB.CONTINUOUS,
                                      name=f'{order} delivery variable', column=None)
                   for order in Order.orders}
payment_variables = {group: mdrp.addVar(lb=0, ub=GRB.INFINITY, obj=0, vtype=GRB.CONTINUOUS,
                                        name=f'{group} payment variable', column=None)
                     for group in Group.groups}
courier_variables = {courier: mdrp.addVar(lb=0, ub=1, obj=0, vtype=GRB.CONTINUOUS,
                                          name=f'{courier} start variable', column=None)
                     for courier in Courier.couriers}

if consider_objective:
    print("Defining objective.")
    if cost_penalty_active:
        mdrp.setObjective(
            quicksum(payment_variables[group] for group in Group.groups)
            + cost_penalty_value * quicksum(1 - order_variables[order] for order in Order.orders)
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
            ),
            name=f'{group} delivery payment constraint')
        for group in Group.groups
    }

    pay_for_time = {
        group: mdrp.addConstr(
            payment_variables[group] >= Data.MIN_PAY_PER_HOUR * group.get_total_on_time() / 60,
            name=f'{group} time payment constraint')
        for group in Group.groups
    }

deliver_all_orders = {
    order: mdrp.addConstr(
        quicksum(fragment_variables[fragment] for fragment in Fragment.fragments_by_order[order])
        == order_variables[order],
        name=f'{order} delivery constraint')
    for order in Order.orders
}

if cost_penalty_active is False:
    all_orders_on = {order: mdrp.addConstr(order_variables[order] == 1,
                                           name=f'{order} activation constraint')
                     for order in Order.orders}
else:
    all_orders_on = {order: mdrp.addConstr(order_variables[order] <= 1,
                                           name=f'{order} activation constraint')
                     for order in Order.orders}

flow_in_equals_flow_out = {
    node: mdrp.addConstr(
        quicksum(
            fragment_variables[fragment] for fragment in Fragment.fragments_by_arrival_node[node]
        )
        == quicksum(
            fragment_variables[fragment]
            for fragment in Fragment.fragments_by_departure_node[node]
        ),
        name=f'{node} flow constraint')
    for node in Node.nodes
    if type(node.location) == Restaurant
}

couriers_start_once = {
    courier: mdrp.addConstr(
        quicksum(
            fragment_variables[fragment]
            for fragment in Fragment.departure_fragments_by_courier[courier]
        )
        == courier_variables[courier],
        name=f'{courier} start constraint')
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
                        >= quicksum(fragment_variables[fragment] for fragment in arc_fragments),
                        name=f'{arc} predecessor VI')
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
                        >= quicksum(fragment_variables[fragment] for fragment in arc_fragments),
                        name=f'{arc} successor VI')
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
    fragment_variables[fragment].vtype = GRB.INTEGER
for order in order_variables:
    order_variables[order].vtype = GRB.BINARY
for courier in courier_variables:
    courier_variables[courier].vtype = GRB.BINARY

mdrp._best_solution_value = GRB.INFINITY
mdrp._solved_subproblems = dict()
mdrp._subproblem_gurobi_time = 0.0


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


def add_lazy_optimality_cut_on_fragments(model: Model, group: Group, fragments: FrozenSet[Fragment],
                                         objective: float) -> None:
    """Add a lazy optimality cut on the given fragments, for the given model."""
    model.cbLazy(
        payment_variables[group] >= objective * (
                quicksum(fragment_variables[fragment] for fragment in fragments) - len(fragments) + 1))


def add_lazy_feasibility_predecessor_cut_on_arcs(model: Model, infeasible_arcs: Set[Arc],
                                                 all_arcs: Set[Arc]) -> None:
    """
    Add a lazy feasibility cut on the given infeasible_arcs, for the given model.

    The sum of the infeasible infeasible_arcs must be less than or equal to the number of
    infeasible infeasible_arcs minus 1, plus the sum of currently non-activated, alternative
    predecessors.
    """
    fragments_for_cut = set(Fragment.get_fragments_from_arcs(list(infeasible_arcs)))
    arc_predecessors = Arc.get_pred_to_arcs(infeasible_arcs)
    alternative_predecessors = set(arc_predecessors) - all_arcs
    alternative_fragments = Fragment.get_fragments_from_arcs(list(alternative_predecessors))
    model.cbLazy(quicksum(fragment_variables[fragment] for fragment in fragments_for_cut)
                 <= len(fragments_for_cut) - 1
                 + quicksum(fragment_variables[fragment] for fragment in alternative_fragments))


def add_lazy_feasibility_cut_on_orders(model: Model, group: Group, feasible_orders: Set[Order],
                                       infeasible_orders: Set[Order]) -> None:
    """Add a lazy feasibility cut on the given orders, for the given model."""
    for inf_order in infeasible_orders:
        orders_for_cut = list(feasible_orders)
        orders_for_cut.append(inf_order)
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
        if suggest_and_repair_solutions:
            calculate_solutions = True
        else:
            calculate_solutions = False

        gurobi_model_objective = model.cbGet(GRB.Callback.MIPSOL_OBJ)

        if gurobi_model_objective > model._best_solution_value + 0.01:
            summary_string = f'{round(gurobi_model_objective, 2)} (rejected)'
            calculate_solutions = False
        else:
            summary_string = f'{round(gurobi_model_objective, 2)} -> '

            # There is some gap between the best found solution and the current
            # solution being investigated. Keep track of this gap, and if the
            # gap decreases below 0, do not suggest this solution to gurobi
            gap_between_found_and_best_solutions = model._best_solution_value - gurobi_model_objective
            if gap_between_found_and_best_solutions < 0.01:
                gap_between_found_and_best_solutions = 0

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
            gurobi_solution_group_objective = model.cbGetSolution(payment_variables[group])
            for i in range(len(group_fragments_solutions)):
                if group_fragments_solutions[i] > 0.9:
                    fragment = group_fragments[i]
                    gurobi_solution_group_fragments.add(fragment)
                    if fragment.departure_location != fragment.arrival_location or len(fragment.order_list) > 0:
                        gurobi_solution_group_arcs.add(group_fragments[i].arc)
                        for order in group_fragments[i].order_list:
                            gurobi_solution_group_orders.add(order)
            if not summary_output:
                print(f'Gurobi assigned orders: {gurobi_solution_group_orders}')

            """
            Check for feasibility and optimality of subproblem.
            
            Put subproblem through the submodel. If comes out infeasible, add
            a feasibility cut. If comes out feasible but under-estimates the
            objective, add an optimality cut. Both cuts are done on arcs.
            
            Once cuts have been made, we can then move on to suggesting and
            repairing solutions. 
            """

            # Solve UFModel on arcs, and add optimality/feasibility cut if necessary
            subproblem = SubModel(group, list(gurobi_solution_group_arcs),
                                  list(gurobi_solution_group_orders))
            subproblem.optimize()
            model._subproblem_gurobi_time += subproblem.getAttr("Runtime")
            if subproblem.getAttr(GRB.Attr.Status) == GRB.OPTIMAL:
                # Add an optimality cut if necessary, else save gurobi's solution
                objective_value = subproblem.getAttr(GRB.Attr.ObjVal)
                if objective_value > gurobi_solution_group_objective + 0.01:
                    # Add an optimality cut
                    if not summary_output:
                        print('Gurobi solution too low, adding optimality cut on arcs')
                    add_lazy_optimality_cut_on_arcs(model, group, frozenset(gurobi_solution_group_arcs),
                                                    objective_value)
                    saved_solution[group] = (objective_value, 0, subproblem.convert_to_timed_path_fragments())
                else:
                    if not summary_output:
                        print('Gurobi solution optimal on arcs')
                    saved_solution[group] = (gurobi_solution_group_objective, 0, gurobi_solution_group_fragments)

            else:
                # subproblem infeasible, add feasibility cut
                assert subproblem.getAttr(GRB.Attr.Status) == GRB.INFEASIBLE
                if not summary_output:
                    print('Gurobi solution infeasible, adding feasibility cut on arcs')
                infeasible_arcs = subproblem._arcs
                # No-good feasibility cut
                model.cbLazy(quicksum(fragment_variables[fragment] for arc in infeasible_arcs for fragment in
                                      Fragment.fragments_by_arc[arc])
                             <= len(infeasible_arcs) - 1)
                saved_solution[group] = (
                    group.get_total_on_time() / 60 * Data.MIN_PAY_PER_HOUR, 10000 * len(gurobi_solution_group_orders),
                    set())

            """
            Suggest master problem solutions to Gurobi.
            
            Only do this if the parameter 'suggest_solutions' is True. If the
            subproblem has already been solved, then we retrieve that solution
            to the subproblem. If it hasn't yet been solved, then we solve it.
            
            Add feasibility and optimality cuts on the orders where necessary.
            """
            if not calculate_solutions:
                continue
            # check if (group, orders) pair already considered
            if (group, frozenset(gurobi_solution_group_orders)) in mdrp._solved_subproblems:
                if not summary_output:
                    print(f'Already solved on {(group, gurobi_solution_group_orders)}, retrieving solution')
                # If group objective value too low, add cut and suggest saved solution
                # If group objective is fine, save existing solution
                subnetwork = model._solved_subproblems[(group, frozenset(gurobi_solution_group_orders))]
                objective = subnetwork.get_best_objective()
                extra_costs = subnetwork.get_undelivered_orders_cost()
                fragments = subnetwork.get_best_fragments()
                saved_solution[group] = (objective, extra_costs, fragments)
                # print("True")

            else:
                # Solve UFModel on orders, and add optimality/feasibility cut if necessary
                if not summary_output:
                    print(f'Looking for better solutions on orders for {group}')
                submodel = SubModel.get_arc_model(group, list(gurobi_solution_group_orders))
                if submodel is None:
                    group_arcs = Arc.get_arcs_with_orders(list(gurobi_solution_group_orders), group)
                    submodel = SubModel(group, group_arcs, list(gurobi_solution_group_orders),
                                        cost_penalty_active=True, time_limit=10, save_model=True,
                                        output_type="None")
                    if not summary_output:
                        print(f'Created new submodel on {gurobi_solution_group_orders}')
                else:
                    if not summary_output:
                        print(f'Retrieved existing submodel on {gurobi_solution_group_orders}')
                    submodel.setParam('TimeLimit', submodel.getParamInfo('TimeLimit')[2] + 10)
                submodel.optimize()

                no_of_solutions = submodel.getAttr(GRB.Attr.SolCount)
                while no_of_solutions == 0:
                    submodel.setParam('TimeLimit', submodel.getParamInfo('TimeLimit')[2] + 10)
                    submodel.optimize()
                    no_of_solutions = submodel.getAttr(GRB.Attr.SolCount)
                model._subproblem_gurobi_time += submodel.getAttr("Runtime")

                if not summary_output:
                    if submodel.getAttr(GRB.Attr.Status) == GRB.OPTIMAL:
                        print('Found optimal solution on group')
                    elif submodel.getAttr(GRB.Attr.Status) == GRB.TIME_LIMIT:
                        print('Time limit reached')
                    else:
                        print(f'Status = {submodel.getAttr(GRB.Attr.Status)}')

                if not summary_output:
                    print(f'Found {no_of_solutions} solutions')
                if no_of_solutions == 0:
                    calculate_solutions = False
                    summary_string += f'rejected (no solutions found for {group})'
                    continue

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
                # print("Hola")
                # print(objective_value, extra_costs, saved_solution[group])
                # print(submodel.get_undelivered_orders())
                # print(submodel._orders)
                if objective_value + extra_costs + 0.01 < saved_solution[group][0] + saved_solution[group][1]:
                    if not summary_output:
                        print(f'Found solution of {objective_value + extra_costs}, suggesting new solution')
                    saved_solution[group] = (objective_value, extra_costs, solution_fragments)
                    # print("Hello")
                    # print(objective_value, extra_costs, solution_fragments)

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

                difference = objective_value + extra_costs - gurobi_solution_group_objective
                if difference > 0.01:
                    gurobi_model_objective += difference
                    # Gap between current solution and best solution decreases
                    gap_between_found_and_best_solutions -= difference
                    if gap_between_found_and_best_solutions < 0.01:
                        calculate_solutions = False
                        summary_string += f'{round(gurobi_model_objective, 2)}+ (rejected)'
                        continue

                # Solved to optimality, save submodel
                if submodel.getAttr(GRB.Attr.Status) == GRB.OPTIMAL:
                    if not summary_output:
                        print(f'Solved submodel to optimality with solution of {objective_value}')
                    mdrp._solved_subproblems[(group, frozenset(gurobi_solution_group_orders))] \
                        = SubNetwork(group, frozenset(gurobi_solution_group_orders), objective_value, extra_costs,
                                     set(solution_fragments))
        if calculate_solutions:
            # Check to see if saved solution better than current best solution
            total_solution_value = 0
            for group in Group.groups:
                total_solution_value += saved_solution[group][0] + saved_solution[group][1]
            for solution_order in Order.orders:
                if model.cbGetSolution(order_variables[solution_order]) < 0.1:
                    total_solution_value += 10000
            summary_string += f'{round(total_solution_value, 2)}'
            if total_solution_value + 0.01 < mdrp._best_solution_value:
                summary_string += ' (saved)'
                if not summary_output:
                    print(f'New solution of {total_solution_value} best found so far')
                mdrp._best_solution_value = total_solution_value
                mdrp._best_solution_fragments: List[Fragment] = list()
                for group in Group.groups:
                    mdrp._best_solution_fragments += saved_solution[group][2]
                mdrp._best_solution_group_values = {group: saved_solution[group][0] for group in Group.groups}
            else:
                if not summary_output:
                    print(f'New solution of {total_solution_value} not best solution')
            if not summary_output:
                print('-' * 50)
        if summary_output:
            print(f't = {get_program_run_time()}, {summary_string}')

    # Suggest saved solution to gurobi
    if where == GRB.Callback.MIPNODE:
        if mdrp._best_solution_value + 0.01 < model.cbGet(GRB.Callback.MIPNODE_OBJBST):
            if not summary_output:
                print(f'\nSuggesting solution of {mdrp._best_solution_value} to Gurobi')
            # print('Best Solution')
            activated_fragments = set(mdrp._best_solution_fragments)
            # print(activated_fragments)
            deactivated_fragments = set(Fragment.fragments) - activated_fragments
            model.cbSetSolution(list(fragment_variables[fragment] for fragment in activated_fragments),
                                list(mdrp._best_solution_fragments.count(fragment) for fragment in activated_fragments))
            model.cbSetSolution(list(fragment_variables[fragment] for fragment in deactivated_fragments),
                                list(0 for _ in range(len(deactivated_fragments))))
            model.cbSetSolution(list(payment_variables[group] for group in mdrp._best_solution_group_values),
                                list(mdrp._best_solution_group_values[group] for group in
                                     mdrp._best_solution_group_values))
            objVal = model.cbUseSolution()
            if not summary_output:
                print(f'Gurobi found objective of {objVal}')
            else:
                print(f'Suggested {round(mdrp._best_solution_value, 2)}')
            if objVal > mdrp._best_solution_value + 0.01:
                print('Infeasible suggested fragment set')
                for fragment in activated_fragments:
                    print(fragment)


# mdrp.setParam('TuneTimeLimit', 36000)
# mdrp.tune()
if group_by_off_time or node_at_order_times or time_discretisation > 1:
    mdrp.optimize(callback)
else:
    mdrp.optimize()
print(f"\nProgram finished running at t = {get_program_run_time()}")
print(f'Gurobi spent {round(mdrp._subproblem_gurobi_time, 2)} seconds in the callback.')

activated_arcs_by_group = {group: list() for group in Group.groups}
orders_by_group = {group: list() for group in Group.groups}
for fragment in fragment_variables:
    if fragment_variables[fragment].x > 0.9:
        if fragment.departure_location != fragment.arrival_location or len(fragment.order_list) > 0:
            group = fragment.group
            activated_arcs_by_group[group].append(fragment.arc)
            for order in fragment.order_list:
                orders_by_group[group].append(order)

itineraries = list()
for group in activated_arcs_by_group:
    submodel = UntimedFragmentsMDRP(group, activated_arcs_by_group[group], orders_by_group[group])
    submodel.optimize()
    itineraries += submodel.get_arc_network_paths()

solution_checker.check_solution(set(Order.orders), itineraries)

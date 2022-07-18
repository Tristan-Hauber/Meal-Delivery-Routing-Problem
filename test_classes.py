"""
Tests relating to the class Location.

Created on Wed Apr 27 14:36:21 2022

@author: Tristan
"""

import unittest
from classes import Data, Location, Courier, Group, Restaurant, Order, Sequence, Arc, Node, Fragment


class ClassesTestCase(unittest.TestCase):
    def setUp(self):
        Location.locations = []
        Courier.couriers = []
        Group.groups = []
        Restaurant.restaurants = []
        Order.orders = []
        Sequence.sequences = set()
        Arc.arcs = set()
        Arc.arcs_by_group_and_order_set_and_arrival_location = dict()
        Node.nodes = []
        Fragment.fragments = []

        self.loc1 = Location(100, 20)
        self.loc2 = Location(54, 60)

        Data.TRAVEL_SPEED = 3
        Data.MAX_CLICK_TO_DOOR = 90
        Data.PICKUP_SERVICE_TIME = 4
        Data.DROPOFF_SERVICE_TIME = 4

        self.cour1 = Courier('c1', 53, 52, 500, 556)
        self.cour2 = Courier('c2', 45, 86, 443, 556)
        self.cour3 = Courier('c3', 39, 50, 245, 875)
        self.cour4 = Courier('c4', 70, 70, 240, 875)

        Group.group_couriers(group_by_off_time=True)
        self.group1 = Group.groups[0]
        self.group2 = Group.groups[1]

        self.rest1 = Restaurant('r1', 40, 40)

        self.ord1 = Order('o1', 40, 60, 450, 'r1', 462)
        self.ord2 = Order('o2', 60, 40, 487, 'r1', 505)
        self.ord3 = Order('o3', 60, 60, 451, 'r1', 478)

        self.seq1 = Sequence([self.ord1], self.rest1)
        self.seq2 = Sequence([self.ord2], self.rest1)
        self.seq3 = Sequence([self.ord3], self.rest1)
        self.seq4 = self.seq1.add_order(self.ord3)
        self.seq5 = Sequence([self.ord3, self.ord1], self.rest1)
        self.seq6 = Sequence([self.ord1, self.ord3, self.ord2], self.rest1)
        self.seq7 = Sequence([self.ord3, self.ord1, self.ord2], self.rest1)
        self.empt_seq = Sequence([], self.rest1)
        self.hom1_seq = Sequence([], self.cour1)
        self.hom2_seq = Sequence([], self.cour2)

        self.arc1 = Arc(self.group1, self.seq1, self.rest1)
        self.arc2 = Arc(self.group1, self.seq2, self.rest1)
        self.arc3 = Arc(self.group1, self.seq3, self.rest1)
        self.arc4 = Arc(self.group1, self.hom2_seq, self.rest1)
        self.arc5 = Arc(self.group1, self.seq6, self.group1)
        self.arc6 = Arc(self.group1, self.seq4, self.group1)
        self.arc7 = Arc(self.group1, self.seq5, self.group1)
        self.arc8 = Arc(self.group1, self.hom1_seq, self.rest1)
        self.wait = Arc(self.group1, self.empt_seq, self.rest1)
        self.arc9 = Arc(self.group2, self.seq1, self.rest1)
        self.arc10 = Arc(self.group1, self.seq4, self.rest1)
        self.arc11 = Arc(self.group1, self.seq2, self.group1)

        self.node1 = Node(self.group1, self.cour1, 500)
        self.node2 = Node(self.group1, self.cour2, 443)
        self.node3 = Node(self.group1, self.rest1, 462)
        self.node4 = Node(self.group1, self.rest1, 478)
        self.node5 = Node(self.group1, self.rest1, 505)
        self.node6 = Node(self.group1, self.group1, 556)
        self.node7 = Node(self.group1, self.rest1, 460)

        self.frag1 = Fragment(self.node1, self.arc8)
        self.frag2 = Fragment(self.node2, self.arc4)
        self.frag3 = Fragment(self.node3, self.wait)
        self.frag4 = Fragment(self.node4, self.arc10)
        self.frag5 = Fragment(self.node5, self.arc11)

    def test_create_location_adds_to_list_of_locations(self):
        self.assertTrue(self.loc1 in Location.locations)

    def test_add_location_has_correct_data(self):
        self.assertEqual(self.loc1.get_location(), (100, 20))

    def test_get_distance_to(self):
        self.assertAlmostEqual(self.loc1.get_distance_to(self.loc2), 60.96, places=1)

    def test_get_time_to(self):
        self.assertEqual(self.loc1.get_time_to(self.loc2), 21)

    def test_create_courier_has_correct_data(self):
        self.assertEqual(self.cour1.code, 'c1')
        self.assertEqual(self.cour1.get_location(), (53, 52))
        self.assertEqual(self.cour1.on_time, 500)
        self.assertEqual(self.cour1.off_time, 556)

    def test_create_courier_adds_to_list(self):
        self.assertTrue(self.cour1 in Courier.couriers)

    def test_group_has_correct_data(self):
        self.assertEqual(self.group1.code, 'g1')
        self.assertEqual(self.group1.off_time, 556)
        self.assertEqual(self.group1.couriers, [self.cour1, self.cour2])

    def test_create_group_adds_to_list(self):
        self.assertTrue(self.group1 in Group.groups)

    def test_get_earliest_arrival_at(self):
        self.assertEqual(self.group1.get_earliest_arrival_at(self.rest1), 461)
        self.assertEqual(self.group2.get_earliest_arrival_at(self.rest1), 251)

    def test_get_total_on_time(self):
        self.assertEqual(self.group1.get_total_on_time(), 169)

    def test_create_restaurant_adds_to_list(self):
        self.assertTrue(self.rest1 in Restaurant.restaurants)

    def test_restaurant_has_correct_data(self):
        self.assertEqual(self.rest1.code, 'r1')
        self.assertEqual(self.rest1.get_location(), (40, 40))

    def test_create_order_adds_to_list(self):
        self.assertTrue(self.ord1 in Order.orders)

    def test_order_has_correct_data(self):
        self.assertEqual(self.ord1.code, 'o1')
        self.assertEqual(self.ord1.get_location(), (40, 60))
        self.assertEqual(self.ord1.placement_time, 450)
        self.assertEqual(self.ord1.departure_restaurant, self.rest1)
        self.assertEqual(self.ord1.earliest_departure_time, 462)
        self.assertEqual(self.ord1.latest_departure_time, 529)

    def test_order_added_to_restaurant_list(self):
        self.assertTrue(self.ord1 in self.rest1.orders)

    # TODO: def test_departure_after_ready(self):

    # TODO: def test_ready_after_placement(self):

    def test_create_sequence_adds_to_list(self):
        self.assertTrue(self.seq1 in Sequence.sequences)

    def test_sequence_has_correct_data(self):
        self.assertEqual(self.seq1.order_list, [self.ord1])
        self.assertEqual(self.seq1.earliest_departure_time, 462)
        self.assertEqual(self.seq1.departure_location, self.rest1)
        self.assertEqual(self.seq1.travel_time, 11)
        self.assertEqual(self.seq1.latest_departure_time, 529)

    def test_add_order_not_modify_original_sequence(self):
        self.assertEqual(self.seq1.order_list, [self.ord1])

    def test_add_order_correct_data(self):
        self.assertEqual(self.seq4.order_list, [self.ord1, self.ord3])
        self.assertEqual(self.seq4.earliest_departure_time, 478)
        self.assertEqual(self.seq4.departure_location, self.rest1)
        self.assertEqual(self.seq4.travel_time, 22)
        self.assertEqual(self.seq4.latest_departure_time, 519)

    def test_sequence_dominates_when_created(self):
        self.assertTrue(self.seq6 in Sequence.sequences)
        self.assertFalse(self.seq7 in Sequence.sequences)

    def test_sequence_does_not_dominate_when_different_ends(self):
        self.assertTrue(self.seq4 in Sequence.sequences)
        self.assertTrue(self.seq5 in Sequence.sequences)

    def test_empty_sequence_has_correct_data(self):
        self.assertEqual(self.empt_seq.order_list, [])
        self.assertEqual(self.empt_seq.earliest_departure_time, 462)
        self.assertEqual(self.empt_seq.departure_location, self.rest1)
        self.assertEqual(self.empt_seq.travel_time, 0)
        self.assertEqual(self.empt_seq.latest_departure_time, 566)

    def test_home_sequence_has_correct_data(self):
        self.assertEqual(self.hom2_seq.order_list, [])
        self.assertEqual(self.hom2_seq.earliest_departure_time, 443)
        self.assertEqual(self.hom2_seq.departure_location, self.cour2)
        self.assertEqual(self.hom2_seq.travel_time, 0)
        self.assertEqual(self.hom2_seq.latest_departure_time, 556)

    # TODO: def test_throws_error_when_add_order_from_wrong_restaurant(self):

    # TODO: def test_throws_error_when_create_sequence_wrong_restaurant(self):

# TODO: def test_throws_error_when_depart_from_group

    def test_create_arc_adds_to_list(self):
        self.assertTrue(self.arc1 in Arc.arcs)

    def test_arc_has_correct_data(self):
        self.assertEqual(self.arc1.group, self.group1)
        self.assertEqual(self.arc1.sequence, self.seq1)
        self.assertEqual(self.arc1.arrival_location, self.rest1)
        self.assertEqual(self.arc1.sequence.order_list, [self.ord1])
        self.assertEqual(self.arc1.departure_location, self.rest1)
        self.assertEqual(self.arc1.earliest_departure_time, 462)
        self.assertEqual(self.arc1.latest_departure_time, 529)
        self.assertEqual(self.arc1.earliest_arrival_time, 484)
        self.assertEqual(self.arc1.latest_arrival_time, 551)
        self.assertEqual(self.arc1.travel_time, 22)

    def test_arc_dominates_when_created(self):
        self.assertTrue(self.arc6 in Arc.arcs)
        self.assertFalse(self.arc7 in Arc.arcs)

    def test_get_pred_to_arc(self):
        pred = self.arc2.get_pred()
        self.assertEqual(pred, {self.arc1, self.arc3, self.arc4})

    def test_get_pred_ignores_arc_with_same_order(self):
        self.assertFalse(self.arc1 in self.arc10.get_pred())

    # TODO: def test_get_pred_ignores_dominated_arcs(self):

    # TODO: def test_get_pred_to_arcs(self):

    # TODO: def test_arr_loc_cannot_be_different_group(self):

    def test_home_arc_has_correct_data(self):
        self.assertEqual(self.arc5.group, self.group1)
        self.assertEqual(self.arc5.sequence, self.seq6)
        self.assertEqual(self.arc5.arrival_location, self.group1)
        self.assertEqual(self.arc5.order_list, [self.ord1, self.ord3, self.ord2])
        self.assertEqual(self.arc5.departure_location, self.rest1)
        self.assertEqual(self.arc5.earliest_departure_time, 505)
        self.assertEqual(self.arc5.latest_departure_time, 519)
        self.assertEqual(self.arc5.earliest_arrival_time, 538)
        self.assertEqual(self.arc5.latest_arrival_time, 552)
        self.assertEqual(self.arc5.travel_time, 33)

    def test_wait_arc_has_correct_data(self):
        self.assertEqual(self.wait.group, self.group1)
        self.assertEqual(self.wait.sequence, self.empt_seq)
        self.assertEqual(self.wait.arrival_location, self.rest1)
        self.assertEqual(self.wait.order_list, [])
        self.assertEqual(self.wait.departure_location, self.rest1)
        self.assertEqual(self.wait.earliest_departure_time, 462)
        self.assertEqual(self.wait.latest_departure_time, 556)
        self.assertEqual(self.wait.earliest_arrival_time, 462)
        self.assertEqual(self.wait.latest_arrival_time, 556)
        self.assertEqual(self.wait.travel_time, 0)

    def test_enter_arc_has_correct_data(self):
        self.assertEqual(self.arc4.group, self.group1)
        self.assertEqual(self.arc4.sequence, self.hom2_seq)
        self.assertEqual(self.arc4.arrival_location, self.rest1)
        self.assertEqual(self.arc4.sequence.order_list, [])
        self.assertEqual(self.arc4.departure_location, self.cour2)
        self.assertEqual(self.arc4.earliest_departure_time, 443)
        self.assertEqual(self.arc4.latest_departure_time, 538)
        self.assertEqual(self.arc4.earliest_arrival_time, 461)
        self.assertEqual(self.arc4.latest_arrival_time, 556)
        self.assertEqual(self.arc4.travel_time, 18)
        self.assertEqual(self.arc8.group, self.group1)
        self.assertEqual(self.arc8.sequence, self.hom1_seq)
        self.assertEqual(self.arc8.arrival_location, self.rest1)
        self.assertEqual(self.arc8.order_list, [])
        self.assertEqual(self.arc8.departure_location, self.cour1)
        self.assertEqual(self.arc8.earliest_departure_time, 500)
        self.assertEqual(self.arc8.latest_departure_time, 548)
        self.assertEqual(self.arc8.earliest_arrival_time, 508)
        self.assertEqual(self.arc8.latest_arrival_time, 556)
        self.assertEqual(self.arc8.travel_time, 8)

    def test_arc_early_departure_constrained_by_arc(self):
        self.assertEqual(self.arc2.earliest_departure_time, 505)

    def test_not_dominate_entry_arcs(self):
        self.assertTrue(self.arc4 in Arc.arcs)
        self.assertTrue(self.arc8 in Arc.arcs)

    def test_add_node(self):
        self.assertTrue(self.node1 in Node.nodes)
        self.assertEqual(self.node3.group, self.group1)
        self.assertEqual(self.node3.location, self.rest1)
        self.assertEqual(self.node3.time, 462)

    def test_add_fragment(self):
        self.assertEqual(self.frag4.departure_node, self.node4)
        self.assertEqual(self.frag4.arc, self.arc10)
        self.assertEqual(self.frag4.arrival_node, self.node5)
        self.assertEqual(self.frag4.group, self.group1)
        self.assertEqual(self.frag4.departure_location, self.rest1)
        self.assertEqual(self.frag4.departure_time, 478)
        self.assertEqual(self.frag4.arrival_location, self.rest1)
        self.assertEqual(self.frag4.arrival_time, 505)

    def test_add_entr_fragment(self):
        self.assertEqual(self.frag1.departure_node, self.node1)
        self.assertEqual(self.frag1.arc, self.arc8)
        self.assertEqual(self.frag1.arrival_node, self.node5)
        self.assertEqual(self.frag1.group, self.group1)
        self.assertEqual(self.frag1.departure_location, self.cour1)
        self.assertEqual(self.frag1.departure_time, 500)
        self.assertEqual(self.frag1.arrival_location, self.rest1)
        self.assertEqual(self.frag1.arrival_time, 505)
        self.assertEqual(self.frag2.departure_node, self.node2)
        self.assertEqual(self.frag2.arc, self.arc4)
        self.assertEqual(self.frag2.arrival_node, self.node7)
        self.assertEqual(self.frag2.group, self.group1)
        self.assertEqual(self.frag2.departure_location, self.cour2)
        self.assertEqual(self.frag2.departure_time, 443)
        self.assertEqual(self.frag2.arrival_location, self.rest1)
        self.assertEqual(self.frag2.arrival_time, 460)

    def test_add_exit_fragment(self):
        self.assertEqual(self.frag5.departure_node, self.node5)
        self.assertEqual(self.frag5.arc, self.arc11)
        self.assertEqual(self.frag5.arrival_node, self.node6)
        self.assertEqual(self.frag5.group, self.group1)
        self.assertEqual(self.frag5.departure_location, self.rest1)
        self.assertEqual(self.frag5.departure_time, 505)
        self.assertEqual(self.frag5.arrival_location, self.group1)
        self.assertEqual(self.frag5.arrival_time, 556)

    def test_add_wait_fragment(self):
        self.assertEqual(self.frag3.departure_node, self.node3)
        self.assertEqual(self.frag3.arc, self.wait)
        self.assertEqual(self.frag3.arrival_node, self.node4)
        self.assertEqual(self.frag3.group, self.group1)
        self.assertEqual(self.frag3.departure_location, self.rest1)
        self.assertEqual(self.frag3.departure_time, 462)
        self.assertEqual(self.frag3.arrival_location, self.rest1)
        self.assertEqual(self.frag3.arrival_time, 478)


if __name__ == '__main__':
    unittest.main()

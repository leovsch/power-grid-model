// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#include <power_grid_model/main_core/input.hpp>
#include <power_grid_model/optimizer/tap_position_optimizer.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <ranges>

namespace pgm_tap = power_grid_model::optimizer::tap_position_optimizer;
namespace main_core = power_grid_model::main_core;

namespace power_grid_model {

namespace {
using TestComponentContainer =
    Container<ExtraRetrievableTypes<Base, Node, Branch, Branch3, Appliance, Regulator>, Line, Link, Node, Transformer,
              ThreeWindingTransformer, TransformerTapRegulator, Source>;
using TestState = main_core::MainModelState<TestComponentContainer>;

TransformerInput get_transformer(ID id, ID from, ID to, BranchSide tap_side) {
    return TransformerInput{.id = id,
                            .from_node = from,
                            .to_node = to,
                            .from_status = 1,
                            .to_status = 1,
                            .u1 = nan,
                            .u2 = nan,
                            .sn = nan,
                            .uk = nan,
                            .pk = nan,
                            .i0 = nan,
                            .p0 = nan,
                            .winding_from = WindingType::wye_n,
                            .winding_to = WindingType::wye_n,
                            .clock = 0,
                            .tap_side = tap_side,
                            .tap_pos = na_IntS,
                            .tap_min = na_IntS,
                            .tap_max = na_IntS,
                            .tap_nom = na_IntS,
                            .tap_size = nan,
                            .uk_min = nan,
                            .uk_max = nan,
                            .pk_min = nan,
                            .pk_max = nan,
                            .r_grounding_from = nan,
                            .x_grounding_from = nan,
                            .r_grounding_to = nan,
                            .x_grounding_to = nan};
}

ThreeWindingTransformerInput get_transformer3w(ID id, ID node_1, ID node_2, ID node_3) {
    return ThreeWindingTransformerInput{
        .id = id,
        .node_1 = node_1,
        .node_2 = node_2,
        .node_3 = node_3,
        .status_1 = 1,
        .status_2 = 1,
        .status_3 = 1,
        .u1 = nan,
        .u2 = nan,
        .u3 = nan,
        .sn_1 = nan,
        .sn_2 = nan,
        .sn_3 = nan,
        .uk_12 = nan,
        .uk_13 = nan,
        .uk_23 = nan,
        .pk_12 = nan,
        .pk_13 = nan,
        .pk_23 = nan,
        .i0 = nan,
        .p0 = nan,
        .winding_1 = WindingType::wye_n,
        .winding_2 = WindingType::wye_n,
        .winding_3 = WindingType::wye_n,
        .clock_12 = 0,
        .clock_13 = 0,
        .tap_side = Branch3Side::side_1,
        .tap_pos = 0,
        .tap_min = 0,
        .tap_max = 0,
        .tap_nom = 0,
        .tap_size = nan,
        .uk_12_min = nan,
        .uk_12_max = nan,
        .uk_13_min = nan,
        .uk_13_max = nan,
        .uk_23_min = nan,
        .uk_23_max = nan,
        .pk_12_min = nan,
        .pk_12_max = nan,
        .pk_13_min = nan,
        .pk_13_max = nan,
        .pk_23_min = nan,
        .pk_23_max = nan,
        .r_grounding_1 = nan,
        .x_grounding_1 = nan,
        .r_grounding_2 = nan,
        .x_grounding_2 = nan,
        .r_grounding_3 = nan,
        .x_grounding_3 = nan,
    };
}

LineInput get_line_input(ID id, ID from, ID to) {
    return LineInput{.id = id,
                     .from_node = from,
                     .to_node = to,
                     .from_status = 1,
                     .to_status = 1,
                     .r1 = nan,
                     .x1 = nan,
                     .c1 = nan,
                     .tan1 = nan,
                     .r0 = nan,
                     .x0 = nan,
                     .c0 = nan,
                     .tan0 = nan,
                     .i_n = nan};
}
TransformerTapRegulatorInput get_regulator(ID id, ID regulated_object, ControlSide control_side) {
    return TransformerTapRegulatorInput{.id = id,
                                        .regulated_object = regulated_object,
                                        .status = 1,
                                        .control_side = control_side,
                                        .u_set = nan,
                                        .u_band = nan,
                                        .line_drop_compensation_r = nan,
                                        .line_drop_compensation_x = nan};
}

} // namespace

TEST_SUITE_BEGIN("Automatic Tap Changer");
TEST_CASE("Test Transformer ranking") {
    // Minimum test grid
    TestState state;
    std::vector<NodeInput> nodes{{0, 150e3}, {1, 10e3}, {2, 10e3}, {3, 10e3}, {4, 10e3},
                                 {5, 50e3},  {6, 10e3}, {7, 10e3}, {8, 10e3}, {9, 10e3}};
    main_core::add_component<Node>(state, nodes.begin(), nodes.end(), 50.0);

    std::vector<TransformerInput> transformers{
        get_transformer(11, 0, 1, BranchSide::from), get_transformer(12, 0, 1, BranchSide::from),
        get_transformer(13, 5, 7, BranchSide::from), get_transformer(14, 2, 3, BranchSide::from),
        get_transformer(15, 8, 9, BranchSide::from)};
    main_core::add_component<Transformer>(state, transformers.begin(), transformers.end(), 50.0);

    std::vector<ThreeWindingTransformerInput> transformers3w{get_transformer3w(16, 0, 4, 5)};
    main_core::add_component<ThreeWindingTransformer>(state, transformers3w.begin(), transformers3w.end(), 50.0);

    std::vector<LineInput> lines{get_line_input(17, 3, 6), get_line_input(18, 3, 9)};
    main_core::add_component<Line>(state, lines.begin(), lines.end(), 50.0);

    std::vector<LinkInput> links{{19, 2, 1, 1, 1}, {20, 6, 4, 1, 1}, {21, 8, 7, 1, 1}};
    main_core::add_component<Link>(state, links.begin(), links.end(), 50.0);

    std::vector<SourceInput> sources{{22, 0, 1, 1.0, 0, nan, nan, nan}};
    main_core::add_component<Source>(state, sources.begin(), sources.end(), 50.0);

    std::vector<TransformerTapRegulatorInput> regulators{
        get_regulator(23, 11, ControlSide::from), get_regulator(24, 12, ControlSide::from),
        get_regulator(25, 13, ControlSide::from), get_regulator(26, 14, ControlSide::from),
        get_regulator(27, 15, ControlSide::from), get_regulator(28, 16, ControlSide::side_1)};
    main_core::add_component<TransformerTapRegulator>(state, regulators.begin(), regulators.end(), 50.0);

    state.components.set_construction_complete();

    // Subcases
    SUBCASE("Building the graph") {
        using pgm_tap::unregulated_idx;

        // reference graph creation
        // Inserted in order of transformer, transformer3w, line and link
        std::vector<std::pair<Idx, Idx>> expected_edges;
        expected_edges.insert(expected_edges.end(), {{0, 1}, {0, 1}, {5, 7}, {2, 3}, {8, 9}});
        expected_edges.insert(expected_edges.end(), {{0, 4}, {4, 5}, {5, 4}, {0, 5}});
        expected_edges.insert(expected_edges.end(), {{3, 6}, {6, 3}, {3, 9}, {9, 3}});
        expected_edges.insert(expected_edges.end(), {{2, 1}, {1, 2}, {6, 4}, {4, 6}, {8, 7}, {7, 8}});

        pgm_tap::TrafoGraphEdgeProperties expected_edges_prop;
        expected_edges_prop.insert(expected_edges_prop.end(),
                                   {{{3, 0}, 1}, {{3, 1}, 1}, {{3, 2}, 1}, {{3, 3}, 1}, {{3, 4}, 1}});
        expected_edges_prop.insert(expected_edges_prop.end(),
                                   {{{4, 0}, 1}, {unregulated_idx, 1}, {unregulated_idx, 1}, {unregulated_idx, 1}});
        expected_edges_prop.insert(expected_edges_prop.end(), 10, {unregulated_idx, 0});

        std::vector<pgm_tap::TrafoGraphVertex> const expected_vertex_props{{true},  {false}, {false}, {false}, {false},
                                                                           {false}, {false}, {false}, {false}, {false}};

        pgm_tap::TransformerGraph actual_graph = pgm_tap::build_transformer_graph(state);
        pgm_tap::TrafoGraphEdgeProperties actual_edges_prop;

        boost::graph_traits<pgm_tap::TransformerGraph>::vertex_iterator vi, vi_end;
        for (boost::tie(vi, vi_end) = vertices(actual_graph); vi != vi_end; ++vi) {
            CHECK(actual_graph[*vi].is_source == expected_vertex_props[*vi].is_source);
        }

        BGL_FORALL_EDGES(e, actual_graph, pgm_tap::TransformerGraph) { actual_edges_prop.push_back(actual_graph[e]); }

        std::sort(actual_edges_prop.begin(), actual_edges_prop.end());
        std::sort(expected_edges_prop.begin(), expected_edges_prop.end());
        CHECK(actual_edges_prop == expected_edges_prop);
    }

    SUBCASE("Automatic tap unsupported tap side at LV") {
        TestState bad_state;
        std::vector<NodeInput> bad_nodes{{0, 50e3}, {1, 10e3}};
        main_core::add_component<Node>(bad_state, bad_nodes.begin(), bad_nodes.end(), 50.0);

        std::vector<TransformerInput> bad_trafo{get_transformer(2, 0, 1, BranchSide::to)};
        main_core::add_component<Transformer>(bad_state, bad_trafo.begin(), bad_trafo.end(), 50.0);

        std::vector<TransformerTapRegulatorInput> bad_regulators{get_regulator(3, 2, ControlSide::from)};

        CHECK_THROWS_AS(main_core::add_component<TransformerTapRegulator>(bad_state, bad_regulators.begin(),
                                                                          bad_regulators.end(), 50.0),
                        AutomaticTapCalculationError);
    }

    SUBCASE("Process edge weights") {
        // Dummy graph
        pgm_tap::TrafoGraphEdges edge_array = {{0, 1}, {0, 2}, {2, 3}};
        pgm_tap::TrafoGraphEdgeProperties edge_prop{{{0, 1}, 1}, {{-1, -1}, 2}, {{2, 3}, 3}};
        std::vector<pgm_tap::TrafoGraphVertex> vertex_props{{true}, {false}, {false}, {false}};

        pgm_tap::TransformerGraph g{boost::edges_are_unsorted_multi_pass, edge_array.cbegin(), edge_array.cend(),
                                    edge_prop.cbegin(), 4};

        // Vertex properties can not be set during graph creation
        boost::graph_traits<pgm_tap::TransformerGraph>::vertex_iterator vi, vi_end;
        for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
            g[*vi].is_source = vertex_props[*vi].is_source;
        }

        pgm_tap::TrafoGraphEdgeProperties const regulated_edge_weights = get_edge_weights(g);
        pgm_tap::TrafoGraphEdgeProperties const ref_regulated_edge_weights{{{0, 1}, 0}, {{2, 3}, 2}};
        CHECK(regulated_edge_weights == ref_regulated_edge_weights);
    }

    SUBCASE("Sorting transformer edges") {
        pgm_tap::TrafoGraphEdgeProperties const trafoList{
            {Idx2D{1, 1}, pgm_tap::infty}, {Idx2D{1, 2}, 5}, {Idx2D{1, 3}, 4}, {Idx2D{2, 1}, 4}};

        pgm_tap::RankedTransformerGroups const referenceList{{Idx2D{1, 3}, Idx2D{2, 1}}, {Idx2D{1, 2}}, {Idx2D{1, 1}}};

        pgm_tap::RankedTransformerGroups const sortedTrafoList = pgm_tap::rank_transformers(trafoList);
        REQUIRE(sortedTrafoList.size() == referenceList.size());
        for (Idx idx = 0; idx < static_cast<Idx>(sortedTrafoList.size()); ++idx) {
            CAPTURE(idx);
            CHECK(sortedTrafoList[idx] == referenceList[idx]);
        }
    }

    SUBCASE("Ranking complete the graph") {
        pgm_tap::RankedTransformerGroups order = pgm_tap::rank_transformers(state);
        pgm_tap::RankedTransformerGroups const ref_order{
            {Idx2D{3, 0}, Idx2D{3, 1}, Idx2D{4, 0}}, {Idx2D{3, 3}, Idx2D{3, 2}}, {Idx2D{3, 4}}};
        CHECK(order == ref_order);
    }
}

TEST_CASE("Test Tap position optimizer" * doctest::skip(true)) {
    // TODO: Implement unit tests for the tap position optimizer
}
TEST_SUITE_END();
} // namespace power_grid_model
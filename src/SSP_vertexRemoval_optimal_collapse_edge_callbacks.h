// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SSP_VERTEXREMOVAL_OPTIMAL_COLLAPSE_EDGE_CALLBACKS_H
#define SSP_VERTEXREMOVAL_OPTIMAL_COLLAPSE_EDGE_CALLBACKS_H
#include <decimate_func_types.h>
#include <igl/quadric_binary_plus_operator.h>

#include <Eigen/Core>
#include <functional>
#include <vector>
#include <tuple>
#include <set>
  // Prepare callbacks for decimating edges using the qslim optimal placement
  // metric.
  //
  // Inputs:
  //   E  #E by 2 list of working edges
  //   quadrics  reference to list of working per vertex quadrics 
  //   v1  working variable to maintain end point of collapsed edge
  //   v2  working variable to maintain end point of collapsed edge
  // Outputs
  //   cost_and_placement  callback for evaluating cost of edge collapse and
  //     determining placement of vertex (see collapse_edge)
  //   pre_collapse  callback before edge collapse (see collapse_edge)
  //   post_collapse  callback after edge collapse (see collapse_edge)
  void SSP_vertexRemoval_optimal_collapse_edge_callbacks(
    Eigen::MatrixXi & E,
    std::vector<std::tuple<Eigen::MatrixXd,Eigen::RowVectorXd,double> > & 
      quadrics,
    int & v1,
    int & v2,
    decimate_cost_and_placement_func & cost_and_placement,
    decimate_pre_collapse_func & pre_collapse,
    decimate_post_collapse_func & post_collapse);
#endif

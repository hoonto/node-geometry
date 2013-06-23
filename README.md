node-geometry
=============

emscripten port of GEOMETRY C library to Javascript for Node and browser environments

Original C version:
http://people.sc.fsu.edu/~jburkardt/c_src/geometry/geometry.html

Note: result of some methods differ in the Node version than in the C version, please use the "runtest" method as shown in test.js to see the output and compare with the original C output contained in geometry_prb_output.txt

Functions
=========

List of Functions:

* angle_box_2d "boxes" an angle defined by three points in 2D.
* angle_contains_ray_2d determines if an angle contains a ray, in 2D.
* angle_deg_2d returns the angle in degrees swept out between two rays in 2D.
* angle_half_2d finds half an angle in 2D.
* angle_rad_2d returns the angle in radians swept out between two rays in 2D.
* angle_rad_3d returns the angle between two vectors in 3D.
* angle_rad_nd returns the angle between two vectors in ND.
* angle_turn_2d computes a turning angle in 2D.
* anglei_deg_2d returns the interior angle in degrees between two rays in 2D.
* anglei_rad_2d returns the interior angle in radians between two rays in 2D.
* annulus_area_2d computes the area of a circular annulus in 2D.
* annulus_sector_area_2d computes the area of an annular sector in 2D.
* annulus_sector_centroid_2d computes the centroid of an annular sector in 2D.
* r8_atan computes the inverse tangent of the ratio Y / X.
* ball_unit_sample_2d picks a random point in the unit ball in 2D.
* ball_unit_sample_3d picks a random point in the unit ball in 3D.
* ball_unit_sample_nd picks a random point in the unit ball in ND.
* basis_map_3d computes the matrix which maps one basis to another.
* box_01_contains_point_2d reports if a point is contained in the unit box in 2D.
* box_01_contains_point_nd reports if a point is contained in the unit box in ND.
* box_contains_point_2d reports if a point is contained in a box in 2D.
* box_contains_point_nd reports if a point is contained in a box in ND.
* box_ray_int_2d: intersection ( box, ray ) in 2D.
* box_segment_clip_2d uses a box to clip a line segment in 2D.
* circle_arc_point_near_2d : nearest point on a circular arc.
* circle_area_2d computes the area of a circle in 2D.
* circle_dia2imp_2d converts a diameter to an implicit circle in 2D.
* circle_exp_contains_point_2d determines if an explicit circle contains a point in 2d.
* circle_exp2imp_2d converts a circle from explicit to implicit form in 2D.
* circle_imp_contains_point_2d determines if an implicit circle contains a point in 2d.
* circle_imp_line_par_int_2d: ( implicit circle, parametric line ) intersection in 2d.
* circle_imp_point_dist_2d: distance ( implicit circle, point ) in 2D.
circle_imp_point_dist_signed_2d: signed distance ( implicit circle, point ) in 2d.
* circle_imp_point_near_2d: nearest ( implicit circle, point ) in 2D.
* circle_imp_points_2d returns N equally spaced points on an implicit circle in 2D.
* circle_imp_points_3d returns points on an implicit circle in 3D.
* circle_imp_points_arc_2d returns N points on an arc of an implicit circle in 2D.
* circle_imp_print_2d prints an implicit circle in 2D.
* circle_imp_print_2d prints an implicit circle in 3D.
* circle_imp2exp_2d converts a circle from implicit to explicit form in 2D.
* circle_llr2imp_2d converts a circle from LLR to implicit form in 2D.
* circle_lune_area_2d returns the area of a circular lune in 2D.
* circle_lune_centroid_2d returns the centroid of a circular lune in 2D.
* circle_ppr2imp_3d converts a circle from PPR to implicit form in 3D.
* circle_ppr2imp_2d converts a circle from PPR to implicit form in 2D.
* circle_sector_area_2d computes the area of a circular sector in 2D.
* circle_sector_centroid_2d returns the centroid of a circular sector in 2D.
* circle_sector_contains_point_2d : is a point inside a circular sector?
* circle_sector_print_2d prints a circular sector in 2D.
* circle_triangle_area_2d returns the area of a circle triangle in 2D.
* circle_triple_angle_2d returns an angle formed by three circles in 2D.
* circles_imp_int_2d: finds the intersection of two implicit circles in 2D.
* cone_area_3d computes the surface area of a right circular cone in 3D.
* cone_centroid_3d returns the centroid of a cone in 3D.
* cone_volume_3d computes the volume of a right circular cone in 3D.
* conv3d converts 3D data to a 2D projection.
* cos_deg returns the cosine of an angle given in degrees.
* cot_deg returns the cotangent of an angle given in degrees.
* cot_rad returns the cotangent of an angle.
* csc_deg returns the cosecant of an angle given in degrees.
* cube_shape_3d describes a cube in 3D.
* cube_size_3d gives "sizes" for a cube in 3D.
* cylinder_point_dist_3d determines the distance from a cylinder to a point in 3D.
* cylinder_point_dist_signed_3d: signed distance from cylinder to point in 3D.
* cylinder_point_inside_3d determines if a cylinder contains a point in 3D.
* cylinder_point_near_3d: nearest point on a cylinder to a point in 3D.
* cylinder_sample_3d samples a cylinder in 3D.
* cylinder_volume_3d determines the volume of a cylinder in 3D.
* degrees_to_radians converts an angle from degrees to radians.
* dge_det computes the determinant of a matrix factored by SGE_FA.
* dge_fa factors a general matrix.
* dge_sl solves a system factored by SGE_FA.
* direction_pert_3d randomly perturbs a direction vector in 3D.
* direction_uniform_2d picks a random direction vector in 2D.
* direction_uniform_3d picks a random direction vector in 3D.
* direction_uniform_nd generates a random direction vector in ND.
* disk_point_dist_3d determines the distance from a disk to a point in 3D.
* dms_to_radians converts an angle from degrees/minutes/seconds to radians.
* dodec_shape_3d describes a dodecahedron in 3D.
* dodec_size_3d gives "sizes" for a dodecahedron in 3D.
* dual_shape_3d constructs the dual of a shape in 3D.
* dual_size_3d determines sizes for a dual of a shape in 3D.
* ellipse_area_2d returns the area of an ellipse in 2D.
* ellipse_point_dist_2d finds the distance from a point to an ellipse in 2D.
* ellipse_point_near_2d finds the nearest point on an ellipse in 2D.
* ellipse_points_2d returns N points on an tilted ellipse in 2D.
* ellipse_points_arc_2d returns N points on a tilted elliptical arc in 2D.
* enorm0_nd computes the Euclidean norm of a (X-Y) in N space.
* get_seed returns a random seed for the random number generator.
* glob2loc_3d converts from a global to a local coordinate system in 3D.
* halfplane_contains_point_2d reports if a half-plane contains a point in 2d.
* halfspace_imp_triangle_int_3d: intersection ( implicit halfspace, triangle ) in 3D.
* halfspace_norm_triangle_int_3d: intersection ( normal halfspace, triangle ) in 3D.
* halfspace_triangle_int_3d: intersection ( halfspace, triangle ) in 3D.
* haversine computes the haversine of an angle.
* helix_shape_3d computes points on a helix in 3D.
* hexagon_area_2d returns the area of a regular hexagon in 2D.
* hexagon_contains_point_2d finds if a point is inside a hexagon in 2D.
* hexagon_shape_2d returns points on the unit regular hexagon in 2D.
* hexagon_unit_area_2d returns the area of a unit regular hexagon in 2D.
* hexagon_vertices_2d returns the vertices of the unit hexagon in 2D.
* i4_dedekind_factor computes a function needed for a Dedekind sum.
* i4_dedekind_sum computes the Dedekind sum of two I4's.
* i4_factorial2 computes the double factorial function N!!
* i4_gcd finds the greatest common divisor of two I4's.
* i4_lcm computes the least common multiple of two I4's.
* i4_max returns the maximum of two I4's.
* i4_min returns the smaller of two I4's.
* i4_modp returns the nonnegative remainder of I4 division.
* i4_sign returns the sign of an I4.
* i4_swap switches two I4's.
* i4_uniform returns a scaled pseudorandom I4.
* i4_wrap forces an I4 to lie between given limits by wrapping.
* i4col_compare compares columns I and J of an I4COL
* i4col_find_item searches an I4COL for a given value.
* i4col_find_pair_wrap wrap searches an I4COL for a pair of items.
* i4col_sort_a ascending sorts the columns of an integer array.
* i4col_sorted_unique_count counts unique elements in an ICOL array.
* i4col_swap swaps two columns of an integer array.
* i4mat_print prints an I4MAT, with an optional title.
* i4mat_print_some prints some of an I4MAT.
* i4mat_transpose_print prints an I4MAT, transposed.
* i4mat_transpose_print_some prints some of an I4MAT, transposed.
* i4row_compare compares two rows of an I4ROW.
* i4row_sort_a ascending sorts the rows of an I4ROW.
* i4row_swap swaps two rows of an I4ROW.
* i4vec_copy copies an I4VEC.
* i4vec_heap_d reorders an I4VEC into a descending heap.
* i4vec_indicator_new sets an I4VEC to the indicator vector.
* i4vec_lcm returns the least common multiple of an I4VEC.
* i4vec_print prints an I4VEC.
* i4vec_product multiplies the entries of an I4VEC.
* i4vec_reverse reverses the elements of an I4VEC.
* i4vec_sort_heap_a ascending sorts an I4VEC using heap sort.
* i4vec_sorted_unique finds unique elements in a sorted I4VEC.
* i4vec_uniform_new returns a scaled pseudorandom I4VEC.
* i4vec_zero zeroes an I4VEC.
* i4vec2_compare compares pairs of integers stored in two vectors.
* i4vec2_sort_a ascending sorts a vector of pairs of integers.
* i4vec2_sorted_unique finds unique elements in a sorted I4VEC2.
* icos_shape describes a icosahedron.
* icos_size gives "sizes" for an icosahedron.
* line_exp_is_degenerate_nd finds if an explicit line is degenerate in ND.
* line_exp_normal_2d computes the unit normal vector to a line in 2D.
* line_exp_perp_2d computes a line perpendicular to a line and through a point.
* line_exp_point_dist_2d: distance ( explicit line, point ) in 2D.
* line_exp_point_dist_3d: distance ( explicit line, point ) in 3D.
* line_exp_point_dist_signed_2d: signed distance ( explicit line, point ) in 2D.
* line_exp_point_near_2d computes the point on an explicit line nearest a point in 2D.
* line_exp_point_near_3d: nearest point on explicit line to point in 3D.
* line_exp2imp_2d converts an explicit line to implicit form in 2D.
* line_exp2par_2d converts a line from explicit to parametric form in 2D.
* line_exp2par_3d converts an explicit line into parametric form in 3D.
* line_imp_is_degenerate_2d finds if an implicit point is degenerate in 2D.
* line_imp_point_dist_2d: distance ( implicit line, point ) in 2D.
* line_imp_point_dist_signed_2d: signed distance ( implicit line, point ) in 2D.
* line_imp2exp_2d converts an implicit line to explicit form in 2D.
* line_imp2par_2d converts an implicit line to parametric form in 2D.
* line_par_point_dist_2d: distance ( parametric line, point ) in 2D.
* line_par_point_dist_3d: distance ( parametric line, point ) in 3D.
* line_par_point_near_2d: nearest point on parametric line to point in 2D.
* line_par_point_dist_3d: distance ( parametric line, point ) in 3D.
* line_par2exp_2d converts a parametric line to explicit form in 2D.
* line_par2exp_2d converts a parametric line to explicit form in 3D.
* line_par2imp_2d converts a parametric line to implicit form in 2D.
* lines_exp_angle_3d finds the angle between two explicit lines in 3D.
* lines_exp_angle_nd returns the angle between two explicit lines in ND.
* lines_exp_dist_3d computes the distance between two explicit lines in 3D.
* lines_exp_dist_3d_2 computes the distance between two explicit lines in 3D.
* lines_exp_equal_2d determines if two explicit lines are equal in 2D.
* lines_exp_int_2d determines where two explicit lines intersect in 2D.
* lines_exp_near_3d computes nearest points on two explicit lines in 3D.
* lines_exp_parallel_2d determines if two lines are parallel in 2D.
* lines_exp_parallel_3d determines if two lines are parallel in 3D.
* lines_imp_angle_2d finds the angle between two implicit lines in 2D.
* lines_imp_dist_2d determines the distance between two implicit lines in 2D.
* lines_imp_int_2d determines where two implicit lines intersect in 2D.
* lines_par_angle_2d finds the angle between two parametric lines in 2D.
* lines_par_angle_3d finds the angle between two parametric lines in 3D.
* lines_par_dist_3d finds the distance between two parametric lines in 3D.
* lines_par_int_2d determines where two parametric lines intersect in 2D.
* loc2glob_3d converts from a local to global coordinate system in 3D.
* lvec_print prints a logical vector.
* minabs finds a local minimum of F(X) = A * abs ( X ) + B.
* minquad finds a local minimum of F(X) = A * X^2 + B * X + C.
* octahedron_shape_3d describes an octahedron in 3D.
* octahedron_size_3d returns size information for an octahedron in 3D.
* parabola_ex finds the extremal point of a parabola determined by three points.
* parabola_ex2 finds the extremal point of a parabola determined by three points.
* parallelogram_area_2d computes the area of a parallelogram in 2D.
* parallelogram_area_3d computes the area of a parallelogram in 3D.
* parallelogram_contains_point_2d determines if a point is inside a parallelogram in 2D.
* parallelogram_contains_point_3d determines if a point is inside a parallelogram in 3D.
* parallelogram_point_dist_3d: distance ( parallelogram, point ) in 3D.
* parallelepiped_contains_point_3d determines if a point is inside a parallelepiped in 3D.
* parallelepiped_point_dist_3d: distance ( parallelepiped, point ) in 3D.
* perm_check checks that a vector represents a permutation.
* perm_inv inverts a permutation "in place".
* plane_exp_grid_3d computes points and lines making up a planar grid in 3D.
* plane_exp_point_dist_3d: distance ( explicit plane, point ) in 3D.
* plane_exp_normal_3d finds the normal to an explicit plane in 3D.
* plane_exp_pro2 produces 2D coordinates of points that lie in a plane, in 3D.
* plane_exp_pro3 projects points orthographically onto a plane, in 3D.
* plane_exp_project_3d projects points through a point onto a plane in 3D.
* plane_exp2imp_3d converts an explicit plane to implicit form in 3D.
* plane_exp2normal_3d converts an explicit plane to normal form in 3D.
* plane_imp_is_degenerate_3d is TRUE if an implicit plane is degenerate.
* plane_imp_line_par_int_3d: intersection ( implicit plane, parametric line ) in 3D.
* plane_imp_point_dist_3d: distance ( point, implicit plane ) in 3D.
* plane_imp_point_dist_signed_3d: signed distance ( implicit plane, point) in 3
* plane_imp_point_near_3d: nearest point on a implicit plane to a point in 3D.
* plane_imp_segment_near_3d: nearest ( implicit plane, line segment ) in 3D
* plane_imp_triangle_int_3d: intersection ( implicit plane, triangle ) in 3D.
* plane_imp_triangle_int_add_3d is a utility for PLANE_IMP_TRIANGLE_INT_3D.
* plane_imp_triangle_near_3d: nearest ( implicit plane, triangle ) in 3D.
* plane_imp2exp_3d converts an implicit plane to explicit form in 3D.
* plane_imp2normal_3d converts an implicit plane to normal form in 3D.
* plane_normal_basis_3d finds two perpendicular vectors in a plane in 3D.
* plane_normal_line_exp_int_3d: intersection of plane and line in 3D.
* plane_normal_qr_to_xyz: QR_TO_XYZ coordinates for a normal form plane.
* plane_normal_tetrahedron_intersect intersects a plane and a tetrahedron.
* plane_normal_triangle_int_3d: intersection ( normal plane, triangle ) in 3D.
* plane_normal_uniform_3d generates a random normal plane in 3D.
* plane_normal_uniform_nd generates a random normal plane in ND.
* plane_normal_xyz_to_qr: XYZ to QR coordinates for a normal form plane.
* plane_normal2exp_3d converts a normal plane to explicit form in 3D.
* plane_normal2imp_3d converts a normal form plane to implicit form in 3D.
* planes_imp_angle_3d: dihedral angle between implicit planes in 3D.
* points_avoid_point_naive_2d finds if a point is "far enough" from a set of points in 2D.
* points_bisect_line_imp_2d finds the implicit line bisecting the line between two points in 2D.
* points_bisect_line_par_2d finds the parametric line bisecting the line between two points in 2D.
* points_centroid_2d computes the discrete centroid of a point set in 2D.
* points_colin_2d estimates the colinearity of 3 points in 2D.
* points_colin_3d estimates the colinearity of 3 points in 3D.
* points_dist_2d finds the distance between two points in 2D.
* points_dist_3d finds the distance between two points in 3D.
* points_dist_nd finds the distance between two points in ND.
* points_hull_2d computes the convex hull of a set of nodes in 2D.
* points_plot plots a pointset.
* points_point_near_naive_2d finds the nearest point to a given point in 2D.
* points_point_near_naive_3d finds the nearest point to a given point in 3D.
* points_point_near_naive_nd finds the nearest point to a given point in ND.
* points_points_near_naive_2d finds the nearest point to given points in 2D.
* points_points_near_naive_3d finds the nearest point to given points in 3D.
* polar_to_xy converts polar coordinates to XY coordinates.
* polygon_1_2d integrates the function 1 over a polygon in 2D.
* polygon_angles_2d computes the interior angles of a polygon in 2D.
* polygon_area_2d computes the area of a polygon in 2D.
* polygon_area_2d_2 computes the area of a polygon in 2D.
* polygon_area_3d computes the area of a polygon in 3D.
* polygon_area_3d_2 computes the area of a polygon in 3D.
* polygon_centroid_2d computes the centroid of a polygon in 2D.
* polygon_centroid_2d_2 computes the centroid of a polygon in 2D.
* polygon_centroid_3d computes the centroid of a polygon in 3D.
* polygon_contains_point_2d finds if a point is inside a simple polygon in 2D.
* polygon_contains_point_2d_2 finds if a point is inside a convex polygon in 2D.
* polygon_diameter_2d computes the diameter of a polygon in 2D.
* polygon_expand_2d expands a polygon in 2D.
* polygon_inrad_data_2d determines polygonal data from its inner radius in 2D.
* polygon_is_convex determines whether a polygon is convex in 2D.
* polygon_lattice_area_2d computes the area of a lattice polygon in 2D.
* polygon_normal_3d computes the normal vector to a polygon in 3D.
* polygon_outrad_data_2d determines polygonal data from its outer radius in 2D.
* polygon_side_data_2d determines polygonal data from its side length in 2D.
* polygon_solid_angle_3d calculates the projected solid angle of a 3D plane polygon.
* polygon_x_2d integrates the function X over a polygon in 2D.
* polygon_y_2d integrates the function Y over a polygon in 2D.
* polygon_xx_2d integrates the function X*X over a polygon in 2D.
* polygon_xy_2d integrates the function X*Y over a polygon in 2D.
* polygon_yy_2d integrates the function Y*Y over a polygon in 2D.
* polyhedron_area_3d computes the surface area of a polyhedron in 3D.
* polyhedron_centroid_3d computes the centroid of a polyhedron in 3D.
* polyhedron_contains_point_3d determines if a point is inside a polyhedron.
* polyhedron_volume_3d computes the volume of a polyhedron in 3D.
* polyhedron_volume_3d_2 computes the volume of a polyhedron in 3D.
* polyline_arclength_nd computes the arclength of points on a polyline in ND.
* polyline_index_point_nd evaluates a polyline at a given arclength in ND.
* polyline_length_nd computes the length of a polyline in ND.
* polyline_points_nd computes equally spaced points on a polyline in ND.
* polyloop_arclength_nd computes the arclength of points on a polyloop in ND.
* polyloop_points_nd computes equally spaced points on a polyloop in ND.
* provec projects a vector from M space into N space.
* pyramid_volume_3d computes the volume of a pyramid with square base in 3D.
* quad_area_2d computes the area of a quadrilateral in 2D.
* quad_area2_2d computes the area of a quadrilateral in 2D.
* quad_area_3d computes the area of a quadrilateral in 3D.
* quad_contains_point_2d finds if a point is inside a convex quadrilateral in 2D.
* quad_convex_random returns a random convex quadrilateral.
* quad_point_dist_2d finds the distance from a point to a quadrilateral in 2D.
* quad_point_dist_signed_2d: signed distanct ( quadrilateral, point ) in 2D.
* quad_point_near_2d computes the nearest point on a quadrilateral in 2D.
* quat_conj conjugates a quaternion.
* quat_inv inverts a quaternion.
* quat_mul multiplies two quaternions.
* quat_norm computes the norm of a quaternion.
* r4_abs returns the absolute value of an R4.
* r4_nint returns the nearest integer to an R4.
* r8_abs returns the absolute value of an R8.
* r8_acos computes the arc cosine function, with argument truncation.
* r8_asin computes the arc sine function, with argument truncation.
* r8_epsilon returns the R8 round off unit.
* r8_huge returns a "huge" R8.
* r8_max returns the maximum of two R8s.
* r8_min returns the minimum of two R8's.
* r8_modp returns the nonnegative remainder of R8 division.
* r8_nint returns the nearest integer to an R8.
* r8_normal_01 returns a unit pseudonormal R8.
* r8_pi returns the value of pi.
* r8_sign returns the "sign" of an R8.
* r8_sign_opposite_strict is TRUE if two R8's are strictly of opposite sign.
* r8_swap switches two R8s.
* r8_uniform returns a scaled pseudorandom R8.
* r8_uniform_01 returns a unit pseudorandom R8.
* r82vec_part_quick_a reorders an R82VEC as part of a quick sort.
* r82vec_permute permutes an R82VEC in place.
* r82vec_print prints an R82VEC.
* r82vec_sort_heap_index_a does an indexed heap ascending sort of an R82VEC.
* r82vec_sort_quick_a ascending sorts an R82VEC using quick sort.
* r8mat_copy copies one R8MAT to another.
* r8mat_det_2d computes the determinant of a 2 by 2 R8MAT.
* r8mat_det_3d computes the determinant of a 3 by 3 R8MAT.
* r8mat_det_4d computes the determinant of a 4 by 4 R8MAT.
* r8mat_det_5d computes the determinant of a 5 by 5 R8MAT.
* r8mat_inverse_2d inverts a 2 by 2 R8MAT using Cramer's rule.
* r8mat_inverse_3d inverts a 3 by 3 R8MAT using Cramer's rule.
* r8mat_mv multiplies a matrix times a vector.
* r8mat_print prints an R8MAT, with an optional title.
* r8mat_print_some prints some of an R8MAT.
* r8mat_solve uses Gauss-Jordan elimination to solve an N by N linear system.
* r8mat_solve_2d solves a 2 by 2 linear system using Cramer's rule.
* r8mat_transpose_print prints an R8MAT, transposed.
* r8mat_transpose_print_some prints some of an R8MAT, transposed.
* r8mat_uniform_new returns a scaled pseudorandom R8MAT.
* r8mat_uniform_01 returns a unit pseudorandom R8MAT.
* r8mat_uniform_01_new returns a unit pseudorandom R8MAT.
* r8vec_angle_3d computes the angle between two vectors in 3D.
* r8vec_any_normal returns some normal vector to V1.
* r8vec_bracket searches a sorted R8VEC for successive brackets of a value.
* r8vec_copy copies an R8VEC.
* r8vec_cross_product_2d finds the cross product of a pair of R8VEC's in 2D.
* r8vec_cross_product_affine_2d finds the affine cross product in 2D.
* r8vec_cross_product_3d computes the cross product of two R8VEC's in 3D.
* r8vec_cross_product_affine_3d computes the affine cross product in 3D.
* r8vec_distance returns the Euclidean distance between two R8VEC's.
* r8vec_dot_product computes the dot product of a pair of R8VEC's.
* r8vec_dot_product_affine computes the affine dot product.
* r8vec_eq is true two R8VEC's are equal.
* r8vec_gt == ( A1 > A2 ) for R8VEC's.
* r8vec_lt == ( A1 < A2 ) for R8VEC's.
* r8vec_max returns the value of the maximum element in an R8VEC.
* r8vec_mean returns the mean of an R8VEC.
* r8vec_min returns the value of the minimum element in an R8VEC.
* r8vec_negative_strict: all entries of R8VEC are strictly negative.
* r8vec_norm returns the L2 norm of an R8VEC.
* r8vec_norm_affine returns the affine L2 norm of an R8VEC.
* r8vec_normal_01_new returns a unit pseudonormal R8VEC.
* r8vec_normsq returns the squared L2 norm of an R8VEC.
* r8vec_normsq_affine returns the squared affine L2 norm of an R8VEC.
* r8vec_positive_strict: all entries of R8VEC are strictly positive.
* r8vec_print prints an R8VEC.
* r8vec_print_2d prints a 2D vector.
* r8vec_print_3d prints a 3D vector.
* r8vec_scalar_triple_product finds the scalar triple product in 3D.
* r8vec_swap swaps the entries of two R8VEC's.
* r8vec_uniform_new returns a scaled pseudorandom R8VEC.
* r8vec_uniform_01_new returns a unit pseudorandom R8VEC.
* r8vec_uniform_unit_new generates a random direction vector in ND.
* r8vec_variance returns the variance of an R8VEC.
* r8vec_zero zeroes an R8VEC.
* radec_distance_3d - angular distance, astronomical units, sphere in 3D.
* radec_to_xyz converts right ascension/declination to (X,Y,Z) coordinates.
* radians_to_degrees converts an angle from radians to degrees.
* radians_to_dms converts an angle from radians to degrees/minutes/seconds.
* random_initialize initializes the RANDOM random number generator.
* rotation_axis_vector_3d rotates a vector around an axis vector in 3D.
* rotation_axis2mat_3d converts a rotation from axis to matrix format in 3D.
* rotation_axis2quat_3d converts a rotation from axis to quaternion format in 3D.
* rotation_mat_vector applies a marix rotation to a vector in 3d.
* rotation_mat2axis_3d converts a rotation from matrix to axis format in 3D.
* rotation_mat2quat_3d converts a rotation from matrix to quaternion format in 3D.
* rotation_quat_vector applies a quaternion rotation to a vector in 3d.
* rotation_quat2axis_3d converts a rotation from quaternion to axis format in 3D.
* rotation_quat2mat_3d converts a rotation from quaternion to matrix format in 3D.
* rtp_to_xyz converts (R,Theta,Phi) to (X,Y,Z) coordinates.
* s_len_trim returns the length of a string to the last nonblank.
* sec_deg returns the secant of an angle given in degrees.
* segment_contains_point_1d reports if a line segment contains a point in 1D.
* segment_contains_point_2d reports if a line segment contains a point in 2D.
* segment_point_coords_2d: coordinates of a point on a line segment in 2D.
* segment_point_coords_3d: coordinates of a point on a line segment in 3D.
* segment_point_dist_2d: distance ( line segment, point ) in 2D.
* segment_point_dist_3d: distance ( line segment, point ) in 3D.
* segment_point_near_2d finds the point on a line segment nearest a point in 2D.
* segment_point_near_3d finds the point on a line segment nearest a point in 3D.
* segments_curvature_2d computes the curvature of two line segments in 2D.
* segments_dist_2d computes the distance between two line segments in 2D.
* segments_dist_3d computes the distance between two line segments in 3D.
* segments_dist_3d_old computes the distance between two line segments in 3D.
* segments_int_1d computes the intersection of two line segments in 1D.
* segments_int_2d computes the intersection of two line segments in 2D.
* shape_point_dist_2d: distance ( regular shape, point ) in 2D.
* shape_point_near_2d: nearest point ( regular shape, point ) in 2D.
* shape_print_3d prints information about a polyhedron in 3D.
* shape_ray_int_2d: intersection ( regular shape, ray ) in 2D.
* simplex_lattice_layer_point_next: next simplex lattice layer point.
* simplex_lattice_point_next returns the next simplex lattice point.
* simplex_unit_lattice_point_nd: count lattice points.
* simplex_unit_volume_nd computes the volume of the unit simplex in ND.
* simplex_volume_nd computes the volume of a simplex in ND.
* sin_deg returns the sine of an angle given in degrees.
* sin_power_int evaluates the sine power integral.
* soccer_shape_3d describes a truncated icosahedron in 3D.
* soccer_size_3d gives "sizes" for a truncated icosahedron in 3D.
* sort_heap_external externally sorts a list of items into ascending order.
* sphere_cap_area_2d computes the surface area of a spherical cap in 2D.
* sphere_cap_area_3d computes the surface area of a spherical cap in 3D.
* sphere_cap_area_nd computes the area of a spherical cap in ND.
* sphere_cap_volume_2d computes the volume of a spherical cap in 2D.
* sphere_cap_volume_3d computes the volume of a spherical cap in 3D.
* sphere_cap_volume_nd computes the volume of a spherical cap in ND.
* sphere_dia2imp_3d converts a diameter to an implicit sphere in 3D.
* sphere_distance_xyz computes great circle distances on a sphere.
* sphere_distance1 computes great circle distances on a sphere.
* sphere_distance2 computes great circle distances on a sphere.
* sphere_distance3 computes great circle distances on a sphere.
* sphere_exp_contains_point_3d determines if an explicit sphere contains a point in 3D.
* sphere_exp_point_near_3d finds the nearest point on an explicit sphere to a point in 3D.
* sphere_exp2imp_3d converts a sphere from explicit to implicit form in 3D.
* sphere_exp2imp_nd finds an N-dimensional sphere through N+1 points.
* sphere_imp_area_3d computes the surface area of an implicit sphere in 3D.
* sphere_imp_area_nd computes the surface area of an implicit sphere in ND.
* sphere_imp_contains_point_3d determines if an implicit sphere contains a point in 3D.
* sphere_imp_grid_icos_size sizes an icosahedral grid on a sphere.
* sphere_imp_gridfaces_3d produces a grid of triangles on an implicit sphere in 3D.
* sphere_imp_line_project_3d projects a line onto an implicit sphere in 3D.
* sphere_imp_local2xyz_3d converts local to XYZ coordinates on an implicit sphere in 3D.
* sphere_imp_point_near_3d finds the nearest point on an implicit sphere to a point in 3D.
* sphere_imp_point_project_3d projects a point onto an implicit sphere, in 3D.
* sphere_imp_volume_3d computes the volume of an implicit sphere in 3D.
* sphere_imp_volume_nd computes the volume of an implicit sphere in ND.
* sphere_imp_zone_area_3d computes the surface area of a spherical zone in 3D.
* sphere_imp_zone_volume_3d computes the volume of a spherical zone in 3D.
* sphere_imp2exp_3d converts a sphere from implicit to explicit form in 3D.
* sphere_k computes a factor useful for spherical computations.
* sphere_triangle_angles_to_area computes the area of a spherical triangle.
* sphere_triangle_contains_point determines if a spherical triangle contains a point.
* sphere_triangle_sides_to_angles computes spherical triangle angles.
* sphere_triangle_vertices_to_angles computes the angles of a spherical triangle.
* sphere_triangle_vertices_to_area computes the area of a spherical triangle.
* sphere_triangle_vertices_to_centroid gets a spherical triangle centroid.
* sphere_triangle_vertices_to_orientation seeks the orientation of a spherical triangle.
* sphere_triangle_vertices_to_sides_3d computes spherical triangle sides.
* sphere_unit_area_nd computes the surface area of a unit sphere in ND.
* sphere_unit_area_values returns some areas of the unit sphere in ND.
* sphere_unit_sample_2d picks a random point on the unit sphere (circle) in 2D.
* sphere_unit_sample_3d picks a random point on the unit sphere in 3D.
* sphere_unit_sample_3d_2 is a BAD method for sampling the unit sphere in 3D.
* sphere_unit_sample_nd picks a random point on the unit sphere in ND.
* sphere_unit_sample_nd_2 picks a random point on the unit sphere in ND.
* sphere_unit_sample_nd_3 picks a random point on the unit sphere in ND.
* sphere_unit_volume_nd computes the volume of a unit sphere in ND.
* sphere_unit_volume_values returns some volumes of the unit sphere in ND.
* sphere01_distance_xyz computes great circle distances on a unit sphere.
* sphere01_polygon_area returns the area of a spherical polygon.
* sphere01_polygon_area_karney returns the area of a spherical polygon.
* sphere01_triangle_angles_to_area: area of a spherical triangle on the unit sphere.
* sphere01_triangle_sides_to_angles: angles of spherical triangle on unit sphere.
* sphere01_triangle_vertices_to_angles: angles of spherical triangle on unit sphere.
* sphere01_triangle_vertices_to_area: area of a spherical triangle on unit sphere.
* sphere01_triangle_vertices_to_centroid: centroid of spherical triangle on unit sphere.
* sphere01_triangle_vertices_to_midpoints gets the midsides of a spherical triangle.
* sphere01_triangle_vertices_to_sides_3d: sides of spherical triangle on unit sphere.
* string_2d groups line segments into connected lines in 2D.
* super_ellipse_points_2d returns N points on a tilted superellipse in 2D.
* tan_deg returns the tangent of an angle given in degrees.
* tetrahedron_barycentric_3d returns the barycentric coordinates of a point in 3D.
* tetrahedron_centroid_3d computes the centroid of a tetrahedron in 3D.
* tetrahedron_circumsphere_3d computes the circumsphere of a tetrahedron in 3D.
* tetrahedron_contains_point_3d: a tetrahedron contains a point in 3D.
* tetrahedron_dihedral_angles_3d computes dihedral angles of a tetrahedron.
* tetrahedron_edge_length_3d returns edge lengths of a tetrahedron in 3D.
* tetrahedron_face_angles_3d returns the 12 face angles of a tetrahedron 3D.
* tetrahedron_face_areas_3d returns the 4 face areas of a tetrahedron 3D.
* tetrahedron_insphere_3d finds the insphere of a tetrahedron in 3D.
* tetrahedron_lattice_layer_point_next: next tetrahedron lattice layer point.
* tetrahedron_lattice_point_next returns the next tetrahedron lattice point.
* tetrahedron_quality1_3d: "quality" of a tetrahedron in 3D.
* tetrahedron_quality2_3d: "quality" of a tetrahedron in 3D.
* tetrahedron_quality3_3d computes the mean ratio of a tetrahedron.
* tetrahedron_quality4_3d computes the minimum solid angle of a tetrahedron.
* tetrahedron_rhombic_shape_3d describes a rhombic tetrahedron in 3D.
* tetrahedron_rhombic_size_3d gives "sizes" for a rhombic tetrahedron in 3D.
* tetrahedron_sample_3d returns random points in a tetrahedron.
* tetrahedron_shape_3d describes a tetrahedron in 3D.
* tetrahedron_size_3d gives "sizes" for a tetrahedron in 3D.
* tetrahedron_solid_angles_3d computes solid angles of a tetrahedron.
* tetrahedron_unit_lattice_point_num_3d: count lattice points.
* tetrahedron_volume_3d computes the volume of a tetrahedron in 3D.
* timestamp prints the current YMDHMS date as a time stamp.
* tmat_init initializes the geometric transformation matrix.
* tmat_mxm multiplies two geometric transformation matrices.
* tmat_mxp multiplies a geometric transformation matrix times a point.
* tmat_mxp2 multiplies a geometric transformation matrix times N points.
* tmat_mxv multiplies a geometric transformation matrix times a vector.
* tmat_rot_axis applies an axis rotation to the geometric transformation matrix.
* tmat_rot_vector applies a rotation about a vector to the geometric transformation matrix.
* tmat_scale applies a scaling to the geometric transformation matrix.
* tmat_shear applies a shear to the geometric transformation matrix.
* tmat_trans applies a translation to the geometric transformation matrix.
* torus_area_3d returns the area of a torus in 3D.
* torus_volume_3d computes the volume of a torus in 3D.
* tp_to_xyz converts unit spherical TP coordinates to XYZ coordinates.
* triangle_angles_2d computes the angles of a triangle in 2D.
* triangle_angles_2d_new computes the angles of a triangle in 2D.
* triangle_angles_3d computes the angles of a triangle in 3D.
* triangle_angles_3d_new computes the angles of a triangle in 3D.
* triangle_area_2d computes the area of a triangle in 2D.
* triangle_area_3d computes the area of a triangle in 3D.
* triangle_area_3d_2 computes the area of a triangle in 3D.
* triangle_area_3d_3 computes the area of a triangle in 3D.
* triangle_area_heron computes the area of a triangle using Heron's formula.
* triangle_area_vector_3d computes the area vector of a triangle in 3D.
* triangle_barycentric_2d finds the barycentric coordinates of a point in 2D.
* triangle_centroid_2d computes the centroid of a triangle in 2D.
* triangle_centroid_3d computes the centroid of a triangle in 3D.
* triangle_circumcenter_2d computes the circumcenter of a triangle in 2D.
* triangle_circumcenter_2d_2 computes the circumcenter of a triangle in 2D.
* triangle_circumcenter computes the circumcenter of a triangle in ND.
* triangle_circumcircle_2d computes the circumcircle of a triangle in 2D.
* triangle_circumcircle_2d_2 computes the circumcircle of a triangle in 2D.
* triangle_circumradius_2d computes the circumradius of a triangle in 2D.
* triangle_contains_line_exp_3d finds if a line is inside a triangle in 3D.
* triangle_contains_line_par_3d: finds if a line is inside a triangle in 3D.
* triangle_contains_point_2d_1 finds if a point is inside a triangle in 2D.
* triangle_contains_point_2d_2 finds if a point is inside a triangle in 2D.
* triangle_contains_point_2d_3 finds if a point is inside a triangle in 2D.
* triangle_diameter_2d computes the diameter of a triangle in 2D.
* triangle_edge_length_2d returns edge lengths of a triangle in 2D.
* triangle_gridpoints_2d computes gridpoints within a triangle in 2D.
* triangle_incenter_2d computes the incenter of a triangle in 2D.
* triangle_incircle_2d computes the inscribed circle of a triangle in 2D.
* triangle_inradius_2d computes the inradius of a triangle in 2D.
* triangle_is_degenerate_nd finds if a triangle is degenerate in ND.
* triangle_lattice_layer_point_next: next triangle lattice layer point.
* triangle_lattice_point_next returns the next triangle lattice point.
* triangle_line_imp_int_2d finds where an implicit line intersects a triangle in 2D.
* triangle_orientation_2d determines the orientation of a triangle in 2D.
* triangle_orthocenter_2d computes the orthocenter of a triangle in 2D.
* triangle_point_dist_2d: distance ( triangle, point ) in 2D.
* triangle_point_dist_3d: distance ( triangle, point ) in 3D.
* triangle_point_dist_signed_2d: signed distance ( triangle, point ) in 2D.
* triangle_point_near_2d computes the nearest triangle point to a point in 2D.
* triangle_quality_2d: "quality" of a triangle in 2D.
* triangle_right_lattice_point_num_2d: count lattice points.
* triangle_sample returns random points in a triangle.
* triangle_unit_lattice_point_num_2d: count lattice points.
* triangle_xsi_to_xy_2d converts from barycentric to XY coordinates in 2D.
* triangle_xy_to_xsi_2d converts from XY to barycentric in 2D.
* truncated_octahedron_shape_3d describes a truncated octahedron in 3D.
* truncated_octahedron_size_3d gives "sizes" for a truncated octahedron in 3D.
* tube_2d constructs a "tube" of given width around a path in 2D.
* tuple_next2 computes the next element of an integer tuple space.
* vector_directions_nd returns the direction angles of a vector in ND.
* vector_rotate_2d rotates a vector around the origin in 2D.
* vector_rotate_3d rotates a vector around an axis vector in 3D.
* vector_rotate_base_2d rotates a vector around a base point in 2D.
* vector_separation_2d finds the angular separation between vectors in 2D.
* vector_separation_3d finds the angular separation between vectors in 3D.
* vector_separation_nd finds the angular separation between vectors in ND.
* vector_unit_nd normalizes a vector in ND.
* voxels_dist_l1_3d computes the L1 distance between voxels in 3D.
* voxels_dist_l1_nd computes the L1 distance between voxels in ND.
* voxels_line_3d computes voxels along a line in 3D.
* voxels_region_3d arranges a set of voxels into contiguous regions in 3D.
* voxels_step_3d computes voxels along a line from a given point in 3D.
* xy_to_polar converts XY coordinates to polar coordinates.
* xyz_to_radec converts (X,Y,Z) to right ascension/declination coordinates.
* xyz_to_rtp converts (X,Y,Z) to (R,Theta,Phi) coordinates.
* xyz_to_tp converts (X,Y,Z) to (Theta,Phi) coordinates.




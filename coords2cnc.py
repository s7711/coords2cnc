import math

# ----- Geometry helper functions -----

def line_intersection(p1, p2, p3, p4):
    """
    Compute the intersection point of the lines defined by (p1-p2) and (p3-p4).
    Returns None if the lines are parallel.
    """
    x1, y1 = p1; x2, y2 = p2
    x3, y3 = p3; x4, y4 = p4
    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(denom) < 1e-8:
        return None
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
    return (x1 + t * (x2 - x1), y1 + t * (y2 - y1))

def translate_points(points, translation):
    """
    Translates a list of 2D points by the given translation vector.

    Args:
        points (list of tuple): List of (x, y) points to be translated.
        translation (tuple): The (xt, yt) translation vector.

    Returns:
        list of tuple: New list of translated (x, y) points.
    """
    xt, yt = translation

    translated_points = [(x + xt, y + yt) for x, y in points]

    return translated_points

def rotate_points(points, angle, centre):
    """
    Rotates a list of 2D points by the given angle around the specified centre.

    Args:
        points (list of tuple): List of (x, y) points to be rotated.
        angle (float): Rotation angle in degrees.
        centre (tuple): The (xc, yc) point to rotate around.

    Returns:
        list of tuple: New list of rotated (x, y) points.
    """
    # Convert angle from degrees to radians
    radians = math.radians(angle)
    # Precompute cosine and sine of the angle
    cos_theta = math.cos(radians)
    sin_theta = math.sin(radians)
    xc, yc = centre

    rotated_points = []
    for x, y in points:
        # Translate point to origin
        x_translated = x - xc
        y_translated = y - yc
        # Rotate point
        x_rotated = x_translated * cos_theta - y_translated * sin_theta
        y_rotated = x_translated * sin_theta + y_translated * cos_theta
        # Translate point back
        x_final = x_rotated + xc
        y_final = y_rotated + yc
        rotated_points.append((x_final, y_final))

    return rotated_points

def offset_polygon(polygon, offset):
    """
    Returns an offset (expanded) version of the given polygon.
    The offset is applied outward.
    
    For a counter-clockwise polygon, shifting using the left-hand normal
    will expand the polygon.
    """
    offset_poly = []
    n = len(polygon)
    for i in range(n):
        A = polygon[(i - 1) % n]
        B = polygon[i]
        C = polygon[(i + 1) % n]

        # Compute unit normal for edge A->B.
        dx1, dy1 = B[0] - A[0], B[1] - A[1]
        len1 = math.hypot(dx1, dy1)
        if len1 == 0:
            continue  # Skip degenerate segment
        normal1 = (-dy1 / len1, dx1 / len1)

        # Compute unit normal for edge B->C.
        dx2, dy2 = C[0] - B[0], C[1] - B[1]
        len2 = math.hypot(dx2, dy2)
        if len2 == 0:
            continue
        normal2 = (-dy2 / len2, dx2 / len2)

        # Shift the edge A->B and B->C by the offset amount.
        line1_start = (A[0] + offset * normal1[0], A[1] + offset * normal1[1])
        line1_end   = (B[0] + offset * normal1[0], B[1] + offset * normal1[1])
        line2_start = (B[0] + offset * normal2[0], B[1] + offset * normal2[1])
        line2_end   = (C[0] + offset * normal2[0], C[1] + offset * normal2[1])

        intersect = line_intersection(line1_start, line1_end, line2_start, line2_end)
        # If no proper intersection is found, use fallback
        if intersect is None:
            intersect = (B[0] + offset * normal1[0], B[1] + offset * normal1[1])
        offset_poly.append(intersect)
    return offset_poly

def polygon_perimeter(polygon):
    """Return the total length of a closed polygon's edges."""
    total = 0.0
    n = len(polygon)
    for i in range(n):
        x0, y0 = polygon[i]
        x1, y1 = polygon[(i+1) % n]
        total += math.hypot(x1 - x0, y1 - y0)
    return total

def get_point_along_polygon(polygon, target_distance):
    """
    Walk along the closed polygon from its first point.
    Returns the (x,y) coordinate at the distance target_distance along the path.
    """
    accumulated = 0.0
    n = len(polygon)
    for i in range(n):
        p_start = polygon[i]
        p_end = polygon[(i+1) % n]
        seg_length = math.hypot(p_end[0] - p_start[0], p_end[1] - p_start[1])
        if accumulated + seg_length >= target_distance:
            f = (target_distance - accumulated) / seg_length
            x = p_start[0] + f * (p_end[0] - p_start[0])
            y = p_start[1] + f * (p_end[1] - p_start[1])
            return (x, y)
        accumulated += seg_length
    return polygon[-1]

# ----- Ramp Cutting Helper -----

def cut_offset_polygon(ncfile, polygon, config, effective_offset, pass_label=""):
    """
    Generates G-code for ramped cutting along an offset polygon.
    
    Parameters:
      ncfile: File-like object for G-code output.
      polygon: List of (x, y) tuples defining the desired part boundary.
      config: Dictionary containing machine configuration:
          - tool_diameter_mm, cutting_speed, material_thickness_mm, 
            max_z_step_mm, z_ramp_deg, safe_height_mm, ramp_start_height_mm.
      effective_offset: The offset to apply. For the rough pass this is:
            tool_offset + xy_final_cut. For the finishing pass it is tool_offset.
      pass_label: A text label (e.g., "Rough Pass" or "Finishing Pass") that will be printed in the G-code comments.
    """
    try:
        adjusted_polygon = offset_polygon(polygon, effective_offset)
        if len(adjusted_polygon) < 3:
            raise ValueError("Offset polygon invalid (fewer than 3 vertices).")
    except Exception as e:
        ncfile.write(f"(WARNING: Polygon offset failed: {str(e)})\n")
        return

    total_depth = config['material_thickness_mm']
    max_step = config['max_z_step_mm']
    passes = math.ceil(total_depth / max_step)
    
    ramp_start = config['ramp_start_height_mm']  # For example, 1 mm above the workpiece.
    angle_rad = math.radians(config['z_ramp_deg'])
    tan_angle = math.tan(angle_rad)
    
    perim = polygon_perimeter(adjusted_polygon)
    
    ncfile.write(f"(--- {pass_label} ---)\n")
    current_z = ramp_start  # Start at the ramp start height for the first pass.
    for p in range(passes):
        current_depth = -min((p+1) * max_step, total_depth)  # Negative value for cutting depth.
        drop_required = current_z - current_depth  # Total vertical drop needed for this pass.
        d_ramp = drop_required / tan_angle  # Distance along path needed to achieve drop.

        # Warn if the polygon perimeter is insufficient to finish the ramp.
        if perim < d_ramp:
            ncfile.write(f"(WARNING: Polygon perimeter ({perim:.2f}mm) is shorter than ramp distance ({d_ramp:.2f}mm))\n")
            
        ncfile.write(f"(Pass {p+1}; target full depth = {current_depth} mm)\n")
        # Move to the beginning of the polygon for this pass:
        start_x, start_y = adjusted_polygon[0]
        ncfile.write(f"G00 X{start_x:.3f} Y{start_y:.3f} (Move to start point)\n")
        
        # Begin the ramp move along the polygon.
        accumulated_distance = 0.0
        ncfile.write(f"G01 Z{current_z:.3f} F{config['cutting_speed']} (Start ramp move)\n")

        n = len(adjusted_polygon)
        for i in range(n):
            p_start = adjusted_polygon[i]
            p_end = adjusted_polygon[(i+1) % n]
            seg_length = math.hypot(p_end[0] - p_start[0], p_end[1] - p_start[1])
            new_total = accumulated_distance + seg_length

            # Case 1: Entire segment is within the ramp region.
            if accumulated_distance < d_ramp:
                if new_total <= d_ramp:
                    target_z = current_z - (new_total * tan_angle)
                    # Do not go below current_depth.
                    if target_z < current_depth:
                        target_z = current_depth
                    ncfile.write(f"G01 X{p_end[0]:.3f} Y{p_end[1]:.3f} Z{target_z:.3f} F{config['cutting_speed']} (Ramp move)\n")
                    accumulated_distance = new_total
                else:
                    # Case 2: This segment crosses the point where full depth is reached.
                    distance_to_ramp_end = d_ramp - accumulated_distance
                    f = distance_to_ramp_end / seg_length
                    intermediate_x = p_start[0] + f * (p_end[0] - p_start[0])
                    intermediate_y = p_start[1] + f * (p_end[1] - p_start[1])
                    # At the end of the ramp, Z must equal current_depth.
                    ncfile.write(f"G01 X{intermediate_x:.3f} Y{intermediate_y:.3f} Z{current_depth:.3f} F{config['cutting_speed']} (End ramp)\n")
                    # Remainder of segment at constant full depth.
                    ncfile.write(f"G01 X{p_end[0]:.3f} Y{p_end[1]:.3f} Z{current_depth:.3f} F{config['cutting_speed']} (Constant depth)\n")
                    accumulated_distance = new_total
            else:
                # Case 3: Already cutting at constant full depth.
                ncfile.write(f"G01 X{p_end[0]:.3f} Y{p_end[1]:.3f} Z{current_depth:.3f} F{config['cutting_speed']} (Constant depth)\n")
                accumulated_distance = new_total

        # Update the current Z position for the next pass.
        current_z = current_depth

    # Handle the remaining ramp after the final pass
    remaining_ramp_distance = d_ramp
    if remaining_ramp_distance > 0:
        ncfile.write(f"(--- Cutting remaining ramp ---)\n")
        accumulated_distance = 0.0
        n = len(adjusted_polygon)

        for i in range(n):
            p_start = adjusted_polygon[i]
            p_end = adjusted_polygon[(i + 1) % n]
            seg_length = math.hypot(p_end[0] - p_start[0], p_end[1] - p_start[1])
            new_total = accumulated_distance + seg_length

            if new_total >= remaining_ramp_distance:
                # The remaining ramp ends within this segment
                distance_to_ramp_end = remaining_ramp_distance - accumulated_distance
                f = distance_to_ramp_end / seg_length
                intermediate_x = p_start[0] + f * (p_end[0] - p_start[0])
                intermediate_y = p_start[1] + f * (p_end[1] - p_start[1])
                ncfile.write(f"G01 X{intermediate_x:.3f} Y{intermediate_y:.3f} Z{current_depth:.3f} F{config['cutting_speed']} (End of remaining ramp)\n")
                break
            else:
                # Entire segment is part of the remaining ramp
                ncfile.write(f"G01 X{p_end[0]:.3f} Y{p_end[1]:.3f} Z{current_depth:.3f} F{config['cutting_speed']} (Remaining ramp segment)\n")
                accumulated_distance = new_total

    # Retract the tool to safe height after all passes and ramp cleanup are complete.
    ncfile.write(f"G00 Z{config['safe_height_mm']} (Retract tool to safe height)\n\n")

# ----- Two-Pass Cutting Routine -----

def cut_two_pass_polygons(ncfile, polygon, config):
    """
    For increased accuracy, first the polygon is cut with an extra offset (xy_final_cut)
    to remove additional material. Then it's recut at the exact size.
    
    Rough Pass:
      effective_offset = (tool diameter/2) + xy_final_cut
      
    Finishing Pass:
      effective_offset = (tool diameter/2)
    """
    # Rough pass with extra material removal.
    rough_offset = config['tool_diameter_mm'] / 2.0 + config['xy_final_cut']
    ncfile.write("(=== Starting Rough Pass (extra material) ===)\n")
    cut_offset_polygon(ncfile, polygon, config, rough_offset, pass_label="Rough Pass")
    
    # Finishing pass at exact dimensions.
    finish_offset = config['tool_diameter_mm'] / 2.0
    ncfile.write("(=== Starting Finishing Pass (exact size) ===)\n")
    cut_offset_polygon(ncfile, polygon, config, finish_offset, pass_label="Finishing Pass")

def drill_hole(ncfile, holes, drill_diameter, config):
    """
    Generates simplified G-code to drill a circular hole by cutting full circles
    with a ramped descent. The finished hole will have the desired 'drill_diameter'.
    
    The tool’s center follows a circular path whose radius is:
    
      Rmotion = (drill_diameter/2) - (tool_diameter/2)
    
    so that the tool’s edge produces the finished size.
    
    The vertical descent per arc segment is computed automatically from:
    
      computed_zstep = arc_step_mm * tan(z_ramp_deg)
    
    Every full-circle arc (G02) outputs only the necessary I offset (the X offset to the hole’s center)
    and the new Z depth.
    
    Config dictionary must contain:
      - 'tool_diameter_mm'      : Cutting tool diameter.
      - 'safe_height_mm'        : Z height for rapid moves.
      - 'ramp_start_height_mm'  : Height at which ramping begins.
      - 'drill_depth_mm'        : Final drill depth.
      - 'z_ramp_deg'            : Ramp angle (in degrees).
      - 'drill_cutting_speed'   : Feed rate for drilling (mm/min).
    """
    # Extract config parameters.
    tool_diameter   = config['tool_diameter_mm']
    safe_height     = config['safe_height_mm']
    ramp_start      = config['ramp_start_height_mm']
    drill_depth     = config['drill_depth_mm']
    z_ramp_deg      = config['z_ramp_deg']
    drill_speed     = config['drill_cutting_speed']
    
    # Compute the computed vertical step from a fixed horizontal arc step.
    arc_step_mm = 1.0  # fixed horizontal arc travel per command (mm)
    z_ramp_rad = math.radians(z_ramp_deg)
    computed_zstep = arc_step_mm * math.tan(z_ramp_rad)
    
    # Total vertical drop (as a positive number).
    total_drop = ramp_start - drill_depth
    num_steps = int(math.ceil(total_drop / computed_zstep))
    
    # Compute the motion circle radius.
    finished_radius = drill_diameter / 2.0
    tool_radius = tool_diameter / 2.0
    # For internal cutting the tool's center must follow a circle _inside_ the finished hole.
    Rmotion = finished_radius - tool_radius
    
    # Minimal I offset for a full circle: The assumed starting point is chosen as the leftmost point,
    # i.e. at theta = π, which puts the tool at (Xc - Rmotion, Yc). With the arc center at (Xc,Yc),
    # the necessary I offset is:
    #    I_value = Xc - (Xc - Rmotion) = Rmotion.
    I_value = Rmotion  # and we assume J = 0 so we omit it.
    
    # Helper to compute a point on a circle.
    def point_on_circle(center, radius, theta):
        return (center[0] + radius * math.cos(theta),
                center[1] + radius * math.sin(theta))
    
    # Process each hole.
    for (Xc, Yc) in holes:
        # Write a header comment.
        ncfile.write(f"(Drilling hole at X{Xc:.4f}, Y{Yc:.4f}, target diameter = {drill_diameter:.4f})\n")
        # Rapid move to safe height.
        ncfile.write(f"G00 Z{safe_height:.4f} (Move to safe height)\n")
        
        # Choose the starting point on the circle.
        # We choose theta = π so that the starting point is (Xc - Rmotion, Yc).
        start_pt = point_on_circle((Xc, Yc), Rmotion, math.pi)
        ncfile.write(f"G00 X{start_pt[0]:.4f} Y{start_pt[1]:.4f} (Rapid move to circle start)\n")
        
        # Plunge to ramp start height.
        ncfile.write(f"G01 Z{ramp_start:.4f} F{drill_speed:.1f} (Plunge to ramp start height)\n")
        
        current_z = ramp_start
        # For each full-circle pass, lower the tool by computed_zstep.
        for i in range(num_steps):
            target_z = current_z - computed_zstep
            if target_z < drill_depth:
                target_z = drill_depth
            # Issue a simplified full-circle G02 command.
            # (Omitting X, Y, and J since the controller assumes the current XY remains.)
            ncfile.write(f"G02 I{I_value:.4f} Z{target_z:.4f} F{drill_speed:.1f} (Full circle, ramp cut)\n")
            current_z = target_z
            if current_z <= drill_depth:
                break
        
        # Final full circle pass at the final drill depth to complete the cut.
        ncfile.write(f"G02 I{I_value:.4f} Z{drill_depth:.4f} F{drill_speed:.1f} (Final full circle at drill depth)\n")
        
        # Finally, move to the safe height.
        ncfile.write(f"G00 Z{safe_height:.4f} (Move to safe height)\n\n")

# ----- Example Usage -----

if __name__ == "__main__":
    # Machine configuration including the new ramp start height and extra final offset.
    machine_config = {
        'tool_diameter_mm': 3.175,       # Cutting tool diameter
        'cutting_speed': 2000,           # mm/min feed rate
        'drill_cutting_speed': 500,      # Speed while drilling
        'material_thickness_mm': 12,     # Total depth to cut
        'max_z_step_mm': 3,              # Maximum depth per pass
        'z_ramp_deg': 10,                # Ramp angle (degrees)
        'safe_height_mm': 5.0,           # Height for safe moves
        'ramp_start_height_mm': 1.0,     # Z height at which ramping starts (above workpiece)
        'xy_final_cut': 0.5,             # Extra offset material for roughing pass
        'drill_max_z_step_mm': 1.0,
        'drill_depth_mm': -12,
    }
    
    # Define the desired part polygon (final dimensions)
    component_55_polygon1 = [
        (0,0), (0,117), (90.5,117), (120.5,151), (132.5,151),
        (162.5,117), (194,117), (194,265), (235,265), (235,0)
    ]
    component_55_polygon2 = [
        (18,44), (66.2,44), (18,89.8)
    ]
    component_55_polygon3 = [
        (26.8,97), (75,52.2), (75,97)
    ]
    component_55_polygon4 = [
        (99,44), (176.3,44), (99,90)
    ]
    component_55_polygon5 = [
        (110.7,97), (188,51), (188,97)
    ]
    component_55_drillholes_4mm = [
        (6,32), (6,49), (6,79), (6,109),
        (26,32), (26,109),
        (55,32), (55,109),
        (87,32), (87,44), (87,71), (87,97), (87,109),
        (113,32), (113,109),
        (126.5,145),
        (142,32), (142,109),
        (171,32), (171,109),
        (200,32), (200,49), (200,79), (200,109), (200,139), (200,169), (200,199), (200,229), (200,259),
        (229,32)
    ]
    component_55_drillholes_16mm = [
        (16,16), (40.5,16), (69.5,16), (98.5,16), (127.5,16), (156.5,16), (185.5,16), (214.5,16)
    ]
    component_55_drillholes_20mm = [
        (217,64), (217,94), (217,124), (217,154), (217,184), (217,214), (217,244)
    ]

    component_55_polygon1_2 = rotate_points(component_55_polygon1, 180, (0,0))
    component_55_polygon1_2 = translate_points(component_55_polygon1_2, (180, 280))
    component_55_polygon2_2 = rotate_points(component_55_polygon2, 180, (0,0))
    component_55_polygon2_2 = translate_points(component_55_polygon2_2, (180, 280))
    component_55_polygon3_2 = rotate_points(component_55_polygon3, 180, (0,0))
    component_55_polygon3_2 = translate_points(component_55_polygon3_2, (180, 280))
    component_55_polygon4_2 = rotate_points(component_55_polygon4, 180, (0,0))
    component_55_polygon4_2 = translate_points(component_55_polygon4_2, (180, 280))
    component_55_polygon5_2 = rotate_points(component_55_polygon5, 180, (0,0))
    component_55_polygon5_2 = translate_points(component_55_polygon5_2, (180, 280))
    component_55_drillholes_4mm_2 = rotate_points(component_55_drillholes_4mm, 180, (0,0))
    component_55_drillholes_4mm_2 = translate_points(component_55_drillholes_4mm_2, (180, 280))
    component_55_drillholes_16mm_2 = rotate_points(component_55_drillholes_16mm, 180, (0,0))
    component_55_drillholes_16mm_2 = translate_points(component_55_drillholes_16mm_2, (180, 280))
    component_55_drillholes_20mm_2 = rotate_points(component_55_drillholes_20mm, 180, (0,0))
    component_55_drillholes_20mm_2 = translate_points(component_55_drillholes_20mm_2, (180, 280))
    
    with open("part.nc", "w") as ncfile:
        # (The setup g-code would be written before calling these functions.)
        cut_two_pass_polygons(ncfile, component_55_polygon1, machine_config)
        cut_two_pass_polygons(ncfile, component_55_polygon2, machine_config)
        cut_two_pass_polygons(ncfile, component_55_polygon3, machine_config)
        cut_two_pass_polygons(ncfile, component_55_polygon4, machine_config)
        cut_two_pass_polygons(ncfile, component_55_polygon5, machine_config)
        drill_hole(ncfile, component_55_drillholes_4mm, 4.0, machine_config)
        drill_hole(ncfile, component_55_drillholes_16mm, 16.0, machine_config)
        drill_hole(ncfile, component_55_drillholes_20mm, 20.0, machine_config)
        cut_two_pass_polygons(ncfile, component_55_polygon1_2, machine_config)
        cut_two_pass_polygons(ncfile, component_55_polygon2_2, machine_config)
        cut_two_pass_polygons(ncfile, component_55_polygon3_2, machine_config)
        cut_two_pass_polygons(ncfile, component_55_polygon4_2, machine_config)
        cut_two_pass_polygons(ncfile, component_55_polygon5_2, machine_config)
        drill_hole(ncfile, component_55_drillholes_4mm_2, 4.0, machine_config)
        drill_hole(ncfile, component_55_drillholes_16mm_2, 16.0, machine_config)
        drill_hole(ncfile, component_55_drillholes_20mm_2, 20.0, machine_config)



import numpy as np

def generate_ring_mesh_with_boundaries(inner_radius, outer_radius, num_radial, num_angular):
    """
    Generates a mesh for a ring defined between `inner_radius` and `outer_radius`
    in the first quadrant, including quadrilateral elements and boundary elements.

    Parameters:
        inner_radius (float): Inner radius of the ring.
        outer_radius (float): Outer radius of the ring.
        num_radial (int): Number of subdivisions in the radial direction.
        num_angular (int): Number of subdivisions in the angular direction.

    Returns:
        vertices (list): List of (x, y) coordinates of the vertices.
        elements (list): List of quadrilateral elements, defined as lists of vertex indices.
        boundaries (list): List of boundary segments, each defined as [flag, type, v1, v2].
    """
    vertices = []
    elements = []
    boundaries = []

    # Radial and angular increments
    radial_coords = np.linspace(inner_radius, outer_radius, num_radial + 1)
    angular_coords = np.linspace(0, np.pi / 2, num_angular + 1)

    # Generate vertices
    for r in radial_coords:
        for theta in angular_coords:
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            vertices.append((x, y))

    # Generate quadrilateral elements
    for i in range(num_radial):
        for j in range(num_angular):
            # Indices of the vertices in the current quadrilateral
            v1 = i * (num_angular + 1) + j
            v2 = v1 + 1
            v3 = (i + 1) * (num_angular + 1) + j
            v4 = v3 + 1

            # Create a quadrilateral element
            elements.append([v1, v2, v4, v3])

    # Generate boundary elements
    # Inner radius (flag 4)
    for j in range(num_angular):
        v1 = j
        v2 = j + 1
        boundaries.append([4, 1, v1, v2])

    # Outer radius (flag 5)
    offset = num_radial * (num_angular + 1)
    for j in range(num_angular):
        v1 = offset + j
        v2 = offset + j + 1
        boundaries.append([5, 1, v1, v2])

    # Vertical edge at x=0 (flag 1)
    for i in range(num_radial):
        v1 = i * (num_angular + 1)
        v2 = (i + 1) * (num_angular + 1)
        boundaries.append([1, 1, v1, v2])

    # Horizontal edge at y=0 (flag 2)
    for i in range(num_radial):
        v1 = i * (num_angular + 1) + num_angular
        v2 = (i + 1) * (num_angular + 1) + num_angular
        boundaries.append([2, 1, v1, v2])

    return vertices, elements, boundaries


def write_mesh_file(vertices, elements, boundaries, filename="ring_mesh_quad_with_boundaries.mesh"):
    """
    Writes the mesh data to a file in MFEM format, with elements listed before vertices and boundaries.

    Parameters:
        vertices (list): List of (x, y) coordinates of the vertices.
        elements (list): List of quadrilateral elements, defined as lists of vertex indices.
        boundaries (list): List of boundary elements, each defined as [flag, type, v1, v2].
        filename (str): Output filename.
    """
    with open(filename, "w") as f:
        # Write header
        f.write("MFEM mesh v1.0\n\n")
        f.write("dimension\n2\n\n")

        # Write elements
        f.write(f"elements\n{len(elements)}\n")
        for element in elements:
            f.write(f"1 3 {' '.join(map(str, element))}\n")

        # Write boundary elements
        f.write(f"\nboundary\n{len(boundaries)}\n")
        for boundary in boundaries:
            flag, btype, v1, v2 = boundary
            f.write(f"{flag} {btype} {v1} {v2}\n")

        # Write vertices
        f.write(f"\nvertices\n{len(vertices)}\n2\n")
        for vertex in vertices:
            f.write(f"{vertex[0]} {vertex[1]}\n")


# Parameters
inner_radius = 0.9
outer_radius = 1.0
num_radial = 10  # Number of subdivisions in the radial direction
num_angular = 40  # Number of subdivisions in the angular direction

# Generate mesh
vertices, elements, boundaries = generate_ring_mesh_with_boundaries(
    inner_radius, outer_radius, num_radial, num_angular
)

# Write mesh to file
home_dir = "/Users/madisonsheridan/Workspace/Laglos/"
output_filename = home_dir + "data/ring.mesh"
write_mesh_file(vertices, elements, boundaries, output_filename)

print(f"Mesh written to {output_filename}")

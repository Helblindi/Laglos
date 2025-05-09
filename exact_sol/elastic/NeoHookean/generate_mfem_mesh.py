import sys

def generate_mfem_mesh(input_file, output_file):
    # Read vertex locations from input file
    with open(input_file, 'r') as f:
        vertices = [float(line.strip()) for line in f.readlines()]
    
    num_vertices = len(vertices)
    num_elements = num_vertices - 1  # One segment per adjacent pair
    num_boundaries = 2  # First and last vertex as boundary markers

    with open(output_file, 'w') as f:
        f.write("MFEM mesh v1.0\n")
        f.write("\ndimension\n1\n")
        
        # Write elements (one segment per consecutive vertex pair)
        f.write("\nelements\n")
        f.write(f"{num_elements}\n")
        for i in range(num_elements):
            f.write(f"1 1 {i} {i+1}\n")  # Segment type (1), element index pairs

        # Write boundary section with explicit point specification
        f.write("\nboundary\n")
        f.write(f"{num_boundaries}\n")
        f.write(f"1 0 0\n")  # First vertex boundary (point)
        f.write(f"1 0 {num_vertices - 1}\n")  # Last vertex boundary (point)

        # Write vertices section
        f.write("\nvertices\n")
        f.write(f"{num_vertices}\n")
        f.write("1\n")  # 1D space
        for v in vertices:
            f.write(f"{v:.16f}\n")  # High precision to avoid small floating-point shifts

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 generate_mfem_mesh.py input.txt output.mesh")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    generate_mfem_mesh(input_file, output_file)

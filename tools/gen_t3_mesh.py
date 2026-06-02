#!/usr/bin/env python3
"""Generate T3 cantilever beam convergence mesh input files.

Each quad in the Q4 mesh is split into 2 triangles (alternating diagonal direction).
T3 element type = 4.
"""
import os, sys

L, H, t = 10.0, 1.0, 0.1
E, nu = 2.0e11, 0.3
F_total = -100000.0
ETYPE = 4  # T3

meshes = [(2, 1), (4, 2), (8, 4), (16, 8), (32, 16)]

script_dir = os.path.dirname(os.path.abspath(__file__))
proj_dir = os.path.dirname(script_dir)

for M, N in meshes:
    nx = M + 1
    ny = N + 1
    nn = nx * ny
    ne = M * N * 2  # each quad = 2 triangles

    lines = [f"T3 Cantilever Convergence Test {M}x{N} Mesh ({M*N*2} triangles)"]
    lines.append(f"{nn}  1  1  1")

    for j in range(ny):
        y = j * H / N
        for i in range(nx):
            x = i * L / M
            bcx = 1 if i == 0 else 0
            bcy = 1 if i == 0 else 0
            bcz = 1
            nid = j * nx + i + 1
            lines.append(f"{nid}  {bcx}  {bcy}  {bcz}  {x:.6f}  {y:.6f}  0.0")

    f_per_node = F_total / ny
    load_parts = [f"1  {ny}"]
    for j in range(ny):
        nid = j * nx + M + 1
        load_parts.append(f"{nid}  2  {f_per_node:.6f}")
    lines.append("  ".join(load_parts))

    lines.append(f"{ETYPE}  {ne}  1")
    lines.append(f"1  {E:.6e}  {nu:.6f}  {t:.6f}")

    eid = 0
    for j in range(N):
        for i in range(M):
            n1 = j * nx + i + 1
            n2 = j * nx + (i + 1) + 1
            n3 = (j + 1) * nx + (i + 1) + 1
            n4 = (j + 1) * nx + i + 1
            # Two triangles per quad, alternating diagonal
            eid += 1
            if (i + j) % 2 == 0:
                lines.append(f"{eid}  {n1}  {n2}  {n3}  1")
                eid += 1
                lines.append(f"{eid}  {n1}  {n3}  {n4}  1")
            else:
                lines.append(f"{eid}  {n1}  {n2}  {n4}  1")
                eid += 1
                lines.append(f"{eid}  {n2}  {n3}  {n4}  1")

    path = os.path.join(proj_dir, "data", f"t3_conv_{M}x{N}.dat")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Generated {path} ({nn} nodes, {ne} triangles)")

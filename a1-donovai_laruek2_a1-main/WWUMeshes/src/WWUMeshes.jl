module WWUMeshes

export read_obj, write_obj
export gen_mesh, est_normals
export cube_mesh, cylinder_mesh, sphere_mesh, estimate_normals, equivilent, equivilent_tri
export OBJTriangle, OBJMesh

using FileIO
using LinearAlgebra

push!(LOAD_PATH, pwd())

include("GfxBase.jl")
using .GfxBase


""" OBJTriangle
A struct that represents a single triangle in a mesh. """
mutable struct OBJTriangle
    positions::Array{Int, 1} # vertex position indices
    uvs::Array{Int, 1} # vertex texture coordinate indices
    normals::Array{Int, 1} # normal vector indices
end

""" OBJMesh
A struct that represents an indexed triangle mesh for reading from
or writing to OBJ format. """
mutable struct OBJMesh
    positions::Array{Vec3, 1} # all vertex positions
    uvs::Array{Vec2, 1} # all texture coordinates
    normals::Array{Vec3, 1} # all vertex normals
    triangles::Array{OBJTriangle, 1} # the OBJTriangles belonging to the mesh
end

""" read_obj(obj_filename)
Read a mesh in OBJ format from file obj_filename."""
function read_obj(obj_filename)
    m = OBJMesh([], [], [], []) # create a mesh
    open(obj_filename) do f
        for (line_number, line) in enumerate(eachline(f))
            if line == "" || line[1] == "#"
                continue # skip comments
            end
            # Read the line and add its contents to the correct field of m:
            tokens = split(strip(line))
            if tokens[1] == "v" # vertex
                push!(m.positions, Vec3([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "vt" # vertex texture
                push!(m.uvs, Vec2([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "vn" # vertex normal
                push!(m.normals, Vec3([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "f"
                # create a OBJTriangle face:
                points = []
                uvs = []
                normals = []
                # handle faces with no texture and/or normals
                for corner in tokens[2:end]
                    indices = split(corner, '/')
                    if length(indices) == 3 # all 3 present, third is a normal
                        push!(normals, parse(Int, indices[3]))
                    end
                    if length(indices) >= 2 && indices[2] != ""
                        # if there are 2 or more and the second isn't blank, it's a texture
                        push!(uvs, parse(Int, indices[2]))
                    end
                    if length(indices) >= 1 # first value is the position
                        push!(points, parse(Int, indices[1]))
                    else # unless it has none, in which case it's not valid
                        error("in line $line_number: face vertex $corner could not be parsed")
                    end
                end
                # create the triangle and add it to the triangles array
                push!(m.triangles, OBJTriangle(points, uvs, normals))
            end
        end
    end
    return m
end

""" write_obj(obj_filename)
Write the given mesh in OBJ format to file obj_filename."""
function write_obj(obj_filename, mesh::OBJMesh)
    open(obj_filename, "w") do f
        # write all positions:
        for v in mesh.positions
            write(f, "v $(v[1]) $(v[2]) $(v[3])\n")
        end

        # write all texture coords:
        for v in mesh.uvs
            write(f, "vt $(v[1]) $(v[2])\n")
        end
        # write all normals:
        for v in mesh.normals
            write(f, "vn $(v[1]) $(v[2]) $(v[3])\n")
        end

        # write all triangles:
        for tri in mesh.triangles
            write(f, "f $(tri_vertex_str(tri))\n")
        end

    end

end

""" tri_vertex_str(triangle)
Return a string with the indices of applicable positions, texture coordinates,
and normals for a given triangle according to the OBJ specification.
In particular, if p, u, and n are position, vertex and normal, each corner
of the triangle is represented as one of the following:
    p       (position only)
    p/u     (position and texture)
    p//n    (position and normal)
    p/u/n   (position, texture, and normal) """
function tri_vertex_str(triangle::OBJTriangle)
    # determine whether textures and normals are present:
    write_uv = length(triangle.uvs) == length(triangle.positions)
    write_normals = length(triangle.normals) == length(triangle.positions)
    corners = []
    for i = 1:3
        output = "$(triangle.positions[i])"
        if write_uv && !write_normals
            output = output * "/$(triangle.uvs[i])" # * does concatenation(!)
        elseif !write_uv && write_normals
            output = output * "//$(triangle.normals[i])"
        elseif write_uv && write_normals
            output = output * "/$(triangle.uvs[i])/$(triangle.normals[i])"
        end
        push!(corners, output)
    end
    join(corners, " ")
end


""" gen_mesh(outfile, geom, divisionsU, divisionsV)
Generate a mesh and save the result in a file with name outfile.
geom may be "cube", "cylinder", or "sphere".
Cylinder requires divisionsU; sphere requires divisionsU and divisionsV. """
function gen_mesh(outfile, geom, divisionsU=0, divisionsV=0, r=0)
    if geom == "cube"
        mesh = cube_mesh()
    elseif geom == "cylinder"
        mesh = cylinder_mesh(divisionsU)
    elseif geom == "sphere"
        mesh = sphere_mesh(divisionsU, divisionsV)
    elseif geom == "torus"
        mesh = torus_mesh(divisionsU, divisionsV, r)
    end
    write_obj(outfile, mesh)
end


""" est_normals(outfile, infile)
Estimate normals of the mesh stored in infile, saving the result in outfile."""
function est_normals(outfile, infile)
    input_mesh = read_obj(infile)
    mesh = estimate_normals(input_mesh)
    write_obj(outfile, mesh)
end

""" find_point_cyl(o, t, l)
Given an 3d origin o, angle t in rads, and length l, find the coordinates of the end of the
vector originating from o
"""
function find_point_cyl(o, t, l)
    x = o[1] + cos(deg2rad(t)) * l
    y = o[2] 
    z = o[3] + sin(deg2rad(t)) * l
    return Vec3(x, y, z)
end

""" find_point(o, t, l)
Given an 2d origin o, angle t in rads, and length l, find the coordinates of the end of the
vector originating from o
"""
function find_point_circ(o, t, l)
    x = o[1] + cos(deg2rad(t)) * l
    y = o[2] + sin(deg2rad(t)) * l
    return Vec2(x, y)
end


""" cube_mesh()
Return a new OBJMesh representing a 2x2x2 cube centered at the origin and
axis-aligned. """
function cube_mesh()
    positions = []
    uvs = []
    normals = []
    triangles = []
    # key to comments:
    # L/R = x = right/left
    # B/T = y = top/bottom
    # C/F = z = close/far
    push!(positions, Vec3( 1, -1, -1)) # 1 RBC
    push!(positions, Vec3( 1, -1,  1)) # 2 RBF
    push!(positions, Vec3(-1, -1,  1)) # 3 LBF
    push!(positions, Vec3(-1, -1, -1)) # 4 LBC
    push!(positions, Vec3( 1,  1, -1)) # 5 RTC
    push!(positions, Vec3( 1,  1,  1)) # 6 RTF
    push!(positions, Vec3(-1,  1,  1)) # 7 LTF
    push!(positions, Vec3(-1,  1, -1)) # 8 LTC

    # texture coordinates:
    push!(uvs, Vec2(1, 1)) # TR
    push!(uvs, Vec2(0, 1)) # TL
    push!(uvs, Vec2(0, 0)) # BL
    push!(uvs, Vec2(1, 0)) # BR

    # normals:
    push!(normals, Vec3( 1, 0, 0)) # R
    push!(normals, Vec3(-1, 0, 0)) # L
    push!(normals, Vec3( 0, 1, 0)) # U
    push!(normals, Vec3( 0,-1, 0)) # D
    push!(normals, Vec3( 0, 0, 1)) # C
    push!(normals, Vec3( 0, 0,-1)) # F

    # 8 faces, 2 triangles each
    push!(triangles, OBJTriangle([1,2,3], [1,2,3], [4,4,4])) # bottom face 1
    push!(triangles, OBJTriangle([1,3,4], [1,3,4], [4,4,4])) # bottom face 2
    push!(triangles, OBJTriangle([1,5,6], [4,1,2], [1,1,1])) # right face 1
    push!(triangles, OBJTriangle([1,6,2], [4,2,3], [1,1,1])) # right face 2
    push!(triangles, OBJTriangle([2,6,7], [4,1,2], [5,5,5])) # far face 1
    push!(triangles, OBJTriangle([2,7,3], [4,2,3], [5,5,5])) # far face 2
    push!(triangles, OBJTriangle([3,7,8], [2,3,4], [2,2,2])) # left face 1
    push!(triangles, OBJTriangle([3,8,4], [2,4,1], [2,2,2])) # left face 2
    push!(triangles, OBJTriangle([4,8,5], [2,3,4], [6,6,6])) # far face 1
    push!(triangles, OBJTriangle([4,5,1], [2,4,1], [6,6,6])) # far face 2
    push!(triangles, OBJTriangle([5,8,7], [1,2,3], [3,3,3])) # top face 1
    push!(triangles, OBJTriangle([5,7,6], [1,3,4], [3,3,3])) # top face 2

    # julia automatically returns the last value in the function:
    OBJMesh(positions, uvs, normals, triangles)

end


""" cylinder_mesh(n)
Return a new OBJMesh object approximation of a cylinder with radius 1 and
height 2, centered at the origin. The logitudinal axis is aligned with y, and
it is tesselated with n divisions arranged radially around the outer surface.
The ends of the cylinder are disc-shaped caps parallel to the xz plane. See the
assignment writeup for a diagram and details.
"""
function cylinder_mesh(divisionsU)
    positions = []
    uvs = []
    normals = []
    triangles = []

    top_center = Vec3(0, 1, 0)
    bottom_center = Vec3(0, -1, 0)
    int_angle = 360 / divisionsU
    current_angle = 0
    current_delta = 1


    #calculate 2 initial points to form 1 wedge
    push!(positions, top_center) #push to positions[1]
    push!(positions, bottom_center) #push to positions[2]
    push!(positions, find_point_cyl(top_center, current_angle, 1)) #push to positions[3]
    push!(positions, find_point_cyl(bottom_center, current_angle, 1)) #push to positions[4]

    #up/down normals
    push!(normals, Vec3(0, 1, 0)) # U normals[1]
    push!(normals, Vec3(0, -1, 0)) # D normals[2]
    #first two points normals (out normals for vertical)
    push!(normals, find_point_cyl(Vec3(0,0,0), current_angle, 1)) #normals[3]
    push!(normals, find_point_cyl(Vec3(0,0,0), current_angle, 1)) #normals[4]

    push!(uvs, Vec2(.75,.75)) #top circle points uvs[1]
    push!(uvs, Vec2(.25,.75)) #bottom circle points uvs[2]
    push!(uvs, find_point_circ(Vec2(.75,.75), current_angle, .25)) #top  0 circle uvs[3]  
    push!(uvs, find_point_circ(Vec2(.25,.75), current_angle, .25)) #bottom 0 circle uvs[4]
    push!(uvs, Vec2(current_delta, .5)) #uvs[6]
    push!(uvs, Vec2(current_delta, 0)) #uvs[5]

    i = 5
    k = 7
    iter = 1
    while iter <= divisionsU

        current_angle += int_angle
        current_delta -= (1/divisionsU)

        push!(positions, find_point_cyl(top_center, current_angle, 1)) #push to positions[i]
        push!(positions, find_point_cyl(bottom_center, current_angle, 1)) #push to positions[i + 1]
        
        
        push!(uvs, find_point_circ(Vec2(.75,.75), current_angle, .25)) #top circle points at k
        push!(uvs, find_point_circ(Vec2(.25,.75), current_angle, .25)) #bottom circle points at k + 1
        push!(uvs, Vec2(current_delta, .5)) # top side points at k + 2
        push!(uvs, Vec2(current_delta, 0)) #bottom side points at k + 3
        

        #out normals
        push!(normals, find_point_cyl(Vec3(0,0,0), current_angle, 1))
        push!(normals, find_point_cyl(Vec3(0,0,0), current_angle, 1))

        #top face tri:
 
        tri1 = OBJTriangle([i - 2, 1, i], [k - 4, 1, k], [1, 1, 1]) 
        #normal is straight up

        #bottom face tri:
        tri2 = OBJTriangle([i + 1, 2, i - 1], [k - 3, 2, k + 1], [2, 2, 2])
        #normal is straight down

        #side tris:
        tri3 = OBJTriangle([i - 2, i, i - 1], [k - 2, k + 2, k - 1], [i - 2, i, i - 1]) #top right side tri
        tri4 = OBJTriangle([i - 1, i, i + 1], [k - 1, k + 2, k + 3], [i - 1, i, i + 1]) #bottom left side tri 

        push!(triangles, tri1)
        push!(triangles, tri2)
        push!(triangles, tri3)
        push!(triangles, tri4)

        i = i + 2
        k = k + 4
        iter = iter + 1

    end
        
    OBJMesh(positions, uvs, normals, triangles)
end

""" torus_param(R1, R2, u, v)
generates x, y, and z coordinates on the paramaterized surface of a torus with major radius R1, minor radius R2
"""
function torus_param(R1, R2, u, v)
    x = (R1 + R2 * cos(deg2rad(v))) * cos(deg2rad(u))
    z = (R1 + R2 * cos(deg2rad(v))) * sin(deg2rad(u))
    y = R2 * sin(deg2rad(v))

    Vec3(x, y, z)
end

""" torus_mesh(n,m,r)
Generates a torus mesh with n longitudinal divisions, n lateral divisions, and a minor
radius r.
"""
function torus_mesh(n, m, r)
    positions = []
    uvs = []
    normals = []
    triangles = []

    int_current_angle = 0
    ext_current_angle = -180
    int_delta_angle = 360 / n
    ext_delta_angle = 360 / m
    first_loop_done = false

    #interior angle loop
    i = 1
    k = 1
    while i <= n + 1
        j = 1

        #calculate u coordinate
        u = int_current_angle / 360

        #find point interior to the torus to calculate the normal from
        int_point = find_point_cyl(Vec3(0,0,0), int_current_angle, 1)

        ext_current_angle = -180
        #exterior angle loop
        while j <= m + 1

            #find the point on the surface of the torus
            point = torus_param(1, r, int_current_angle, ext_current_angle)
            push!(positions, point)

            #calculate normal to point on surface
            norm = point - int_point
            push!(normals, norm)

            #calculate v based on internal loop angle
            v = (ext_current_angle / 360) + 0.5
            uv = Vec2(u, v)
            push!(uvs, uv)

            if first_loop_done && j > 1
                #start constructing tris
                tri1 = OBJTriangle([k - 1, k - m - 1, k], [k - 1, k - m - 1, k], [k - 1, k - m - 1, k])
                tri2 = OBJTriangle([k - m - 2, k - m - 1, k - 1], [k - m - 2, k - m - 1, k - 1], [k - m - 2, k - m - 1, k - 1])
                push!(triangles, tri1)
                push!(triangles, tri2)
            end

            ext_current_angle += ext_delta_angle
            j += 1
            k += 1
        end
        first_loop_done = true

        int_current_angle += int_delta_angle
        i += 1
    end

    OBJMesh(positions, uvs, normals, triangles)
end


""" sphere_mesh(n, m)
Create a Latitude-Longitude-tesselated approximation of a sphere with radius 1
centered at the origin. There are n divisions around the equator and m
divisions from pole to pole along each line of longitude. The North pole is at
(0,1,0), the South pole at (0,-1,0), and points on the Greenwich meridian are
in the x = 0 plane with z > 0. The u texture coordinate depends on longitude,
with u=0 at 180 degrees West and u=1 at 180 degrees East. The v texture
coordinate varies with latitude with v=0 at the South pole and v=1 at the North
pole. Normals should be normal to the ideal sphere surface. See the assignment
for a diagram and further details. """
function sphere_mesh(n, m)
    # Latitude: m, v, j
    # Longitude: n, u, i

    positions = []
    uvs = []
    normals = []
    triangles = []

    lon_int_angle = 360 / n
    lat_int_angle = 180 / m
    lon_current_angle = 0
    

    first_col_done = false

    i = 0
    j = 0
    k = 0
    while i < n + 1        
        lat_current_angle = 0
        j = 0
        while j < m + 1
            k = k + 1

            #find a paramaterized point on the outside of the sphere
            point = sphere_xyz(lat_current_angle, lon_current_angle)

            #calculate u and v coordinates for that point
            u = lon_current_angle / 360
            v = 1 - ((lat_current_angle) / 180)
            uv = Vec2(u,v)

            #normal vector on a sphere centered at 0,0,0 is equivilent to the position of that point
            push!(positions, point)
            push!(normals, point)
            push!(uvs, uv)

            if(first_col_done)
                #start constructing tris

                if (j == 1)
                    #if at top cap
                    tri1 = OBJTriangle([k, k - 1, k - m - 1], [k, k - 1, k - m - 1], [k, k - 1, k - m - 1]) 
                    push!(triangles, tri1) 
                elseif (j == m)
                    #if at bottom cap
                    tri2 = OBJTriangle([k , k - 1, k - m - 2], [k , k - 1, k - m - 2], [k , k - 1, k - m - 2]) 
                    push!(triangles, tri2)
                else
                    #general - non cap tri squares
                    tri3 = OBJTriangle([k - 1, k - m - 1, k], [k - 1, k - m - 1, k], [k - 1, k - m - 1, k])
                    tri4 = OBJTriangle([k - m - 2, k - m - 1, k - 1], [k - m - 2, k - m - 1, k - 1], [k - m - 2, k - m - 1, k - 1])

                    push!(triangles, tri3)
                    push!(triangles, tri4)
                end


            end

            lat_current_angle = lat_current_angle + lat_int_angle
            j = j + 1
        end
        first_col_done = true

        lon_current_angle = lon_current_angle + lon_int_angle
        i = i + 1
    end

    OBJMesh(positions, uvs, normals, triangles)
end

""" sphere_xyz
given a latitude and longitude, return x,y,z for a circle with radius 1
"""
function sphere_xyz(lat, lon)
    z = cos(deg2rad(lon))sin(deg2rad(lat))
    x = sin(deg2rad(lon))sin(deg2rad(lat))
    y = cos(deg2rad(lat))

    return Vec3(x,y,z)
end 

"""
    estimate_normals(mesh::OBJMesh)
Estimates normals for the given mesh. Overwrites any existing normals and returns a new OBJMesh object.
"""
function estimate_normals(mesh::OBJMesh)
    normals = []

    #initialize empty array to use for normal allocation
    for vertex in mesh.positions
        push!(normals, Vec3(0,0,0))
    end
    v = 1
    for tri in mesh.triangles

        #get indicies of points on our local triangle
        indexA = tri.positions[1]
        indexB = tri.positions[2]
        indexC = tri.positions[3]

        #get actual positions of the triangle verticies
        a = mesh.positions[indexA]
        b = mesh.positions[indexB]
        c = mesh.positions[indexC]

        #calculate the normal
        norm = cross((a - b), (a - c))

        #add to normal to the exisiting normals at each point to average their direction
        normals[indexA] += norm
        normals[indexB] += norm
        normals[indexC] += norm

        #append normal indicies to the local triangle
        tri.normals = [indexA, indexB, indexC]
        
    end

    #normalize normal vectors
    v = 1
    for normal in normals
        normals[v] = normalize(normal)
        v += 1
    end

    OBJMesh(mesh.positions, mesh.uvs, normals, mesh.triangles)
end

"""functon to check if two meshes are the same. should give an explanation if verbose is true, eps is margin of error for float point values"""
function equivilent(file1, file2, verbose, eps::Float64)


    mesh1 = read_obj(file1)
    mesh2 = read_obj(file2)

    error1 = false
    error2 = false
    error3 = false
   
    if size(mesh1.triangles, 1) != size(mesh2.triangles, 1)
      error1 = true
    end

    #check if face in mesh1 exists in mesh2
    error2 = check_list(mesh1.triangles, mesh2.triangles, mesh1, mesh2, eps)

    #check if face in mesh1 exists in mesh2
    error3 = check_list(mesh2.triangles, mesh1.triangles, mesh2, mesh1, eps)

    #error reporting to user
    if error1 || error2 || error3
        if verbose
            if error1
                println("The obj files do not have same number of triangles")
            end
            if error2
                println("a face in mesh1 does not have an equivilent face in mesh2")
            end
            if error3
                println("a face mesh2 does not have an equivilent face in mesh1")
            end 
        else
            println("the two meshes are not equivilent")
        end
    else
        println("the two meshes are equivilent")
    end
end

#func to do looping for equivilancey logic: is element from list1 from mesh1 exist in list2 from mesh2 - 
#returns true/false based on if error is found
function check_list(list1, list2, mesh1::OBJMesh, mesh2::OBJMesh, eps)
    missing = true
    i = 0
    for face1 in list1
        missing = true
        for face2 in list2
            i = i + 1
            equal_check = equivilent_tri(face1, face2, mesh1, mesh2, eps)
            if equal_check
                missing = false
            elseif missing && i == length(list2)
                return true
            end
        end 
    end

    return false
end

"""check if triangle a from mesh1 has the same xyz coodinates as triangle b from mesh2"""
function equivilent_tri(tri_a, tri_b, mesh1::OBJMesh, mesh2::OBJMesh, eps)
    uv_exists = false
    normals_exists = false
    uv_case = -1
    norm_case = -1
    #get triangle a' vertecies
    c = mesh1.positions[tri_a.positions[1]]
    d = mesh1.positions[tri_a.positions[2]]
    e = mesh1.positions[tri_a.positions[3]]

    #get triangle b' vertecies
    f = mesh2.positions[tri_b.positions[1]]
    g = mesh2.positions[tri_b.positions[2]]
    h = mesh2.positions[tri_b.positions[3]]

    #get all xyz values for each point and list order is c_xyz, d_xyz, e_xyz, f_xyz, g_xyz, h_xyz)
    pos_values = Any[]
    pos_values = reorder_array_triplets(pos_values, c, d, e, f, g, h)

    #collect texture locations if needed
    #get triangle a' vertecies
    if length(mesh1.uvs) == length(mesh2.uvs) && length(mesh1.uvs) != 0 
        uvs_exists = true

        uv1 = mesh1.uvs[tri_a.uvs[1]]
        uv2 = mesh1.uvs[tri_a.uvs[2]]
        uv3 = mesh1.uvs[tri_a.uvs[3]]

        #get triangle b' vertecies
        uv4 = mesh2.uvs[tri_b.uvs[1]]
        uv5 = mesh2.uvs[tri_b.uvs[2]]
        uv6 = mesh2.uvs[tri_b.uvs[3]]

        uvs_values = Any[]
        uvs_values = reorder_array_pair(uvs_values, uv1, uv2, uv3, uv4, uv5, uv6)
    end
    
    """collect normal locations if needed
    get triangle a' vertecies """
    if length(mesh1.normals) == length(mesh2.normals) && length(mesh1.normals) != 0 
        normals_exists = true
    
        norm1 = mesh1.normals[tri_a.normals[1]]
        norm2 = mesh1.normals[tri_a.normals[2]]
        norm3 = mesh1.normals[tri_a.normals[3]]

        #get triangle b' vertecies
        norm4 = mesh2.normals[tri_b.normals[1]]
        norm5 = mesh2.normals[tri_b.normals[2]]
        norm6 = mesh2.normals[tri_b.normals[3]]

        norm_values = Any[]
        norm_values = reorder_array_triplets(norm_values, norm1, norm2, norm3, norm4, norm5, norm6)
    end 

    """do calculation on diffrence between xyz values
    check three possible coodinate CCW orderings"""
    position_case = compare_xyz_triplets(pos_values, eps)
    
    #if uv/normals exist, do work
    if uv_exists
        uv_case = compare_uv_pairs(uvs_values, eps)
    end
    if normals_exists
        norm_case = compare_xyz_triplets(norm_values, eps)
    end

    #checks error values and returns if equivilent
    if position_case == 0
    return false
    end

    if (uv_exists && uv_case != position_case) || uv_case == 0
        return false
    end

    if (normals_exists && norm_case != position_case) || norm_case == 0
        return false
    end

    return true
end


"""compares xyz of 6 points given in list(order is c_xyz, d_xyz, e_xyz, f_xyz, g_xyz, h_xyz) and returns
which rotation case we have - if no match, return 0"""
function compare_xyz_triplets(xyz, eps)
    #get x value diff - case1
    dif1 = abs(xyz[4] - xyz[13])
    dif2 = abs(xyz[7] - xyz[16])
    dif3 = abs(xyz[1] - xyz[10])

    #get y value diff - case1
    dif4 = abs(xyz[5] - xyz[14])
    dif5 = abs(xyz[8] - xyz[17])
    dif6 = abs(xyz[2] - xyz[11])

    #get z diff - case 1
    dif7 = abs(xyz[6] - xyz[15])
    dif8 = abs(xyz[9] - xyz[18])
    dif9 = abs(xyz[3] - xyz[12])

    #get x value diff - case2
    dif10 = abs(xyz[4] - xyz[16])
    dif11 = abs(xyz[7] - xyz[10])
    dif12 = abs(xyz[1] - xyz[13])

    #get y value diff - case2
    dif13 = abs(xyz[5] - xyz[17])
    dif14 = abs(xyz[8] - xyz[11])
    dif15 = abs(xyz[2] - xyz[14])

    #get z diff - case2
    dif16 = abs(xyz[6] - xyz[18])
    dif17 = abs(xyz[9] - xyz[12])
    dif18 = abs(xyz[3] - xyz[15])

    #get x value diff - case3
    dif19 = abs(xyz[4] - xyz[10])
    dif20 = abs(xyz[7] - xyz[13])
    dif21 = abs(xyz[1] - xyz[16])

    #get y value diff - case3
    dif22 = abs(xyz[5] - xyz[11])
    dif23 = abs(xyz[8] - xyz[14])
    dif24 = abs(xyz[2] - xyz[17])

    #get z diff - case3
    dif25 = abs(xyz[6] - xyz[12])
    dif26 = abs(xyz[9] - xyz[15])
    dif27 = abs(xyz[3] - xyz[18])

    #check if we have a cirtain CCW rotation case or no matches 
    if dif1 <= eps && dif2 <= eps && dif3 <= eps && dif4 <= eps && dif5 <= eps && dif6 <= eps && dif7 <= eps && dif8 <= eps && dif9 <= eps
        return 1
    elseif dif10 <= eps && dif11 <= eps && dif12 <= eps && dif13 <= eps && dif14 <= eps && dif15 <= eps && dif16 <= eps && dif17 <= eps && dif18 <= eps
        return 2
    elseif dif19 <= eps && dif20 <= eps && dif21 <= eps && dif22 <= eps && dif23 <= eps && dif24 <= eps && dif25 <= eps && dif26 <= eps && dif27 <= eps
        return 3
    else
        return 0
    end
end

"""compare uv texture coodinates: compares uv of 6 points given in list(order is c_uv, d_uv, e_uv, f_uv, g_uv, h_uv) """
function compare_uv_pairs(uvs_values, eps)
    #get u value diff - case1
    dif1 = abs(uvs_values[1] - uvs_values[7])
    dif2 = abs(uvs_values[3] - uvs_values[9])
    dif3 = abs(uvs_values[5] - uvs_values[11])

    #get v value diff - case1
    dif4 = abs(uvs_values[2] - uvs_values[8])
    dif5 = abs(uvs_values[4] - uvs_values[10])
    dif6 = abs(uvs_values[6] - uvs_values[12])

    #get u value diff - case2
    dif7 = abs(uvs_values[1] - uvs_values[9])
    dif8 = abs(uvs_values[3] - uvs_values[11])
    dif9 = abs(uvs_values[5] - uvs_values[7])

    #get v value diff - case2
    dif10 = abs(uvs_values[1] - uvs_values[10])
    dif11 = abs(uvs_values[3] - uvs_values[12])
    dif12 = abs(uvs_values[5] - uvs_values[8])

    #get u value diff - case3
    dif13 = abs(uvs_values[1] - uvs_values[11])
    dif14 = abs(uvs_values[3] - uvs_values[7])
    dif15 = abs(uvs_values[5] - uvs_values[9])

    #get v value diff - case3
    dif16 = abs(uvs_values[1] - uvs_values[12])
    dif17 = abs(uvs_values[3] - uvs_values[8])
    dif18 = abs(uvs_values[5] - uvs_values[10])

    #check if we have a cirtain CCW rotation case or no matches
    if dif1 <= eps && dif2 <= eps && dif3 <= eps && dif4 <= eps && dif5 <= eps && dif6 <= eps 
        return 1
    elseif dif7 <= eps && dif8 <= eps && dif9 <= eps && dif10 <= eps && dif11 <= eps && dif12 <= eps
        return 2
    elseif dif13 <= eps && dif14 <= eps && dif15 <= eps && dif16 <= eps && dif17 <= eps && dif18 <= eps
        return 3
    else
        return 0
    end
end

"""function collects triplet values and puts them into array"""
function reorder_array_triplets(array, a, b, c, d, e, f)
    push!(array, a[1])
    push!(array, a[2])
    push!(array, a[3])
    push!(array, b[1])
    push!(array, b[2])
    push!(array, b[3])
    push!(array, c[1])
    push!(array, c[2])
    push!(array, c[3])
    push!(array, d[1])
    push!(array, d[2])
    push!(array, d[3])
    push!(array, e[1])
    push!(array, e[2])
    push!(array, e[3])
    push!(array, f[1])
    push!(array, f[2])
    push!(array, f[3])

    return array
end
"""function collects pair values and puts them into array"""
function reorder_array_pair(array, a, b, c, d, e, f)
    push!(array, a[1])
    push!(array, a[2])
    push!(array, b[1])
    push!(array, b[2])
    push!(array, c[1])
    push!(array, c[2])
    push!(array, d[1])
    push!(array, d[2])
    push!(array, e[1])
    push!(array, e[2])
    push!(array, f[1])
    push!(array, f[2])

    return array
end

end # module WWUMeshes



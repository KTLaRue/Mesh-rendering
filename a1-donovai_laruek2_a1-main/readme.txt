Torus:

The torus extension was relatively simple, as the body of the code is very similar to a sphere. The main difference was the parameterization of the surface, which can be described by the following function:

	""" torus_param(R1, R2, u, v)
	generates x, y, and z coordinates on the parameterized surface of a torus with major radius R1, minor radius R2
	"""
	function torus_param(R1, R2, u, v)
    		x = (R1 + R2 * cos(deg2rad(v))) * cos(deg2rad(u))
    		z = (R1 + R2 * cos(deg2rad(v))) * sin(deg2rad(u))
    		y = R2 * sin(deg2rad(v))

    		Vec3(x, y, z)
	end

This function maps the x,y,z points on the surface of a torus to the major and minor radii and angles. Using these equations, we were able to associate a single normal and uv coordinate pair to each position in the .obj file, meaning we had an equivalent number of positions, normals, and uvs. This made constructing the shape extremely simple, as the index references for each triangle vertex's position, normal, and uv were exactly the same. This meant there were no special cases for areas of the torus, like there were for the caps of the sphere and cylinder.

The generation of a torus can be invoked in the same way as any other shape, through the gen_mesh function. As directed in the instructions, the arguments are longitudinal divisions, lateral divisions, and minor radius. An example call is:

gen_mesh("results/torus", "torus", 32, 16, .5)


Equivalency check

This extension started out relatively straight forward, but got complicated rather quickly. The idea of checking if two separate meshes have the same amount of triangles, then investigating if those triangles have the same xyz coordinates (within an epsilon error range) as well as investigating the uv values and normal vectors if they exist.

There is definitely an easier and straight forward way of doing this since the implementation ended up being 7 methods long but any other implementation was just causing indexing issues.
The general methodology is that the number of triangles need to be equivalent, then every face in mesh1 must exist in mesh2 and every face in mesh2 must exist in mesh1. Technically, this means that for equivalent faces their xyz coordinates, uv coordinates and normal vectors must be the same if they exist.

The methods to do this check are at the bottom of our code and can be evoked from the julia environment with this general command: equivilent("path/to/mesh1", "path/to/mesh2", verbose, eps). Verbose is a boolean to give you more information about any errors that exist and eps is your desired level of error for a vertex value to be considered equal with another.

# stereo-tracer
This is a rendering program written in C++ that uses Ray Tracing.
User defined scenes are placed in the scenes folder, and output bmps as well as a bmp of the depth of the first object the ray hits will be placed there.

Stereoscopic 3D images are rendered by defining the distance the 'eyes' will be offset from each other.
On default, the image will be a Cross view image for easy viewing (left eye's image is on the right, and vice versa)
Using a negative value for the StereoDist value will produce a Parellel View image (left eye's image is on the left).
Using 0.0 for the value will create a normal single-image render.

NOTE that the aspect ratio is for the entire final file, and not for the half images for each eye.

History

Authored by Max Wei and Mikhail Bessmeltsev

Originally a CPSC 314 final assignment, I added Refraction support and then Stereoscopic 3D support.

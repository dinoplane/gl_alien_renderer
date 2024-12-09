# Alien Gl Renderer 

Cloth Simulation Video here!
[![Link to Video](https://www.youtube.com/watch?v=ppvQZppc3sY)]

The cloth simulation source code is [![here](https://github.com/dinoplane/gl_alien_renderer/blob/main/src/cloth_system.cpp)]


This project is a sandbox for the for me to improve and learn numerous graphics concepts:

I post my progress here: https://x.com/isacsac2017

Below is a running list of finished, in progress, and cut features, in addition to my ramblings.

## Finished features
- Loading glTF Models
- Loading glTF Materials 
- GPU Frustum Culling + standalone models
- Debug Frustums on cameras
- Multiple Cameras
- Post-processing
- Instancing
- Indirect Draws using glMultiDrawElementsIndirect
- Multiple Cameras
- ImGui interface
- Tracy Profiler
- Scene loading
- Cloth Simulation

## (Fixed) Bugs:
- Materials not loading:
    - a result of unable to create a "buffer of textures"
    - fix: use a texture atlas

- Instancing not working with MultiDrawIndirect?
    - This was caused when I used a struct to contain both the visibleInstIndices and the visibleInstCount. 
    - Using sizeof the struct yielded 32 bytes only. sizeof did not account for the fact that there was a vector inside of itself. 
    - I just used a single vector like before to hold the count as well (since it's a uint as well). 
    - I may change this solution if there is other data (that has a different type) to be accounted for. 
- cirno's eyes disappear sometimes... figure out why
    - This was also happening with the cactuar
    - Fragments get discarded because a comparison in the fragment shader keeps firing
        - this is because the alpha cuttoff was greater than the alpha
        - this was caused by the model loading a material out of the range of its material vector
            - out of range exception wasn't thrown because [] operator has **undefined behavior** for indices out of range
        - so the alpha cutoff was garbage sometimes and 0 most of the time.
        - Fixed by making sure I assume that there is a default material in every model.
- Indirect Draw Debugging 
    - First symptom: only one model is drawn
        - each instance is probably using the same model matrix, this was because I was using gl_DrawID and DrawElementsIndirect, so glDrawID was always 0
    - Second symptom: meshes with multiple primitives are not drawn correctly.
        - check the draw calls in renderdoc.
            - seems that the baseVertex and firstIndex members in the IndirectDrawCommand struct s were wrong (I flipped them on the CPU side)
        - even after flipping them something was wrong...
            - after checking the values in the node primitives properties struct in renderdoc, it seems that the struct's alignment was wrong. this is because I used alignas(16) in the struct and the struct only had 2 4 byte members
            - remember, the std430 layout has more flexible rules on struct alignment than std140


Capstone TODOs:
- Implement an L system
- Add concepts in animation
    - Mass spring dampers
    - Particle System
    - Splines
    - Inverse kinematics
- Cloth Simulation
    - Try your best ot make it real time
- Particle Systems
- Animation

TODO
- Abstract organize the logic better
- A renderpass system
- lighting
- baked lighting + instancing... damn im gonna regret this aint i

- Get rid of include clutter
- shadows
- reflections
- animation
- Integrating OpenUSD
- multidraws/indirect draws/culling that actually reduces compute
    - a draw a model instance with 1 draw call
        - not possible thats not what multidraw means
        - instead batch the draw calls (every N is a new instance)
        - we can only really do it this way since if i were to have instance count be greater than one in the draw command, i would need to move around amtrices into a new buffer (inefficient...)
        - the draw command probably needs the index of the model though
    - Steps:
        - 1. Modify the renderer to use indirect draws, DrawElementsIndirect.
            - the glviewer example code has a way to set this up per mesh.
            - very intuitive code, every primitive begins with a draw call
            - a mesh contains a buffer of primitives
            - DrawElementsIndirect is called on the location in the buffer (takes the 20 bytes from that pointer)
        - 2. Aggregate all the primitive's vertices and indices into 1 giant buffer
            - make sure i can render indirectly like that
        - 3. Use MultiDrawElementsIndirect instead, this time to draw everything in that buffer
        - 4. Modify the compute shader so that it modifies the draw out cmd buffer
        - 5. Modify the shaders so that it uses materials and textures (using an index to query them from ssbos);

- deferred shading
- add doxygen
- 
- [Craig Reynold's Boids](https://www.cs.toronto.edu/~dt/siggraph97-course/cwr87/)
- Realistic Ocean Waves with help from [GPU Gems](https://developer.nvidia.com/gpugems/gpugems/part-i-natural-effects/chapter-1-effective-water-simulation-physical-models) and this [video](https://www.youtube.com/watch?v=PH9q0HNBjT4)

- Improve the waves. I made it so that the light can be seen underwater, but I'm not exactly sure if it's correct.
- A new boid behavior? Not sure what I should do.
- A debug interface. I should play around with compute shaders so that I can see the normals of the objects I render. Not sure how to implement the window though.
- Sounds (they do improve the experience...)


Bugs:
- fix noninstanced model shader (im lazy so i dwanna right now)



CUT: 
- Make sure to delete all buffers after exiting // eeeeeeehhh?
    There's no point if memory is not in demand in the lifetime.
    The program will release all its resources upon exit anyway!
    My excuse for allowing memory leaks lol!
- Use git submodules? Debatable
    Don't wanna mess with this, in production code, using public modules is tied to legality issues, 



But yea, I hope I can continue playing around with this!

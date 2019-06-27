require 'strict'

-- define some convenience functions

function get_sphere_inertia (mass, radius)
    local val
    val = 2. / 5. * mass * radius * radius
    return {
        {val, 0.0, 0.0},
        {0.0, val, 0.0},
        {0.0, 0.0, val}
    }
end

function get_box_inertia (mass, width, height, depth)
    local valx, valy, valz
    valx = 1. / 12. * mass * (height*height + depth*depth)
    valy = 1. / 12. * mass * (width*width   + height*height)
    valz = 1. / 12. * mass * (width*width   + depth*depth)
    return {
        {valx,  0.0,  0.0},
        { 0.0, valy,  0.0},
        { 0.0,  0.0, valz}
    }
end

-- define some constants

cart_w = 0.5 -- x-axis
cart_d = 0.2 -- y-axis
cart_h = 0.2 -- z-axis
cart_m = 10.
cart_com = {
    0.0, 0.0, -cart_h/4.0
}

pend_l = 0.5
pend_r = 0.1
pend_m = 1.0
pend_com = {
    0.0, 0.0, pend_l
}

-- build model

bodies = {
    cart     = {mass = cart_m, com = cart_com, inertia = get_box_inertia (cart_m, cart_w, cart_h, cart_d)},
    pendulum = {mass = pend_m, com = pend_com, inertia = get_sphere_inertia (pend_m, pend_r)},
}

meshes = {
    cart = {
        name = "cart",
        dimensions = {cart_w, cart_d, cart_h},
        color = {0.0, 1.0, 0.0},
        mesh_center = {0, 0, 0},
        translate = {0, 0, 0},
        rotate = {0.0, 0.0, 0.0},
        src = "unit_cube.obj",
    },
    cart_com = {
        name = "cart_com",
        dimensions = {0.05, 0.05, 0.05},
        color = {0.0, 0.0, 0.0},
        mesh_center = cart_com,
        translate = {0, 0, 0},
        rotate = {0.0, 0.0, 0.0},
        src = "unit_cube.obj",
    },
    link = {
        name = "link",
        dimensions = {0.02, 0.02, pend_l},
        color = {0.5, 0.5, 0.0},
        mesh_center = {0, 0, 0},
        translate = {0, 0, pend_l/2.0},
        rotate = {0.0, 0.0, 0.0},
        src = "unit_cube.obj",
    },
    pendulum = {
        name = "pendulum",
        dimensions = {pend_r, pend_r, pend_r},
        color = {1.0, 0.0, 0.0},
        mesh_center = {0, 0, 0},
        translate = {0, 0, pend_l},
        rotate = {0.0, 0.0, 0.0},
        src = "unit_sphere_medres.obj",
    },
    pend_com = {
        name = "pend_com",
        dimensions = {0.05, 0.05, 0.05},
        color = {0.0, 0.0, 0.0},
        mesh_center = pend_com,
        translate = {0, 0, 0},
        rotate = {0.0, 0.0, 0.0},
        src = "unit_cube.obj",
    },
}

joints = {
    fixed = {},
    trans_x = {
        { 0., 0., 0., 1., 0., 0.},
    },
    rot_y = {
        { 0., 1., 0., 0., 0., 0.},
    },
}

model = {
    gravity = { 0., 0., -9.81},

    configuration = {
        axis_front = { 1.,  0.,  0.},
        axis_right = { 0., -1.,  0.},
        axis_up    = { 0.,  0.,  1.}
    },

    frames = {
        -- acctuated cart
        {
            name = "cart",
            parent = "ROOT",
            body = bodies.cart,
            joint = joints.trans_x,
            joint_frame = {
                r = {0.5, 0.5, 0.7},
                E = {
                    {1.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0},
                    {0.0, 0.0, 1.0},
                }
            },
            visuals = {
                meshes.cart,
                meshes.cart_com,
            }
        },
        -- unactuated pendulum
        {
            name = "pendulum",
            parent = "cart",
            body = bodies.pendulum,
            joint = joints.rot_y,
            joint_frame = {
                r = {0.0, 0.0, 0.0},
                E = {
                    {1.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0},
                    {0.0, 0.0, 1.0},
                }
            },
            visuals = {
                meshes.pendulum,
                meshes.pend_com,
                meshes.link,
            }
        },
        -- dummy frame for curve plotting
        {
            name = "pendulum_com",
            parent = "pendulum",
            joint = joints.fixed,
            joint_frame = {
                r = pend_com,
                E = {
                    {1.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0},
                    {0.0, 0.0, 1.0},
                }
            },
            visuals = {
            }
        },
    }
}

return model

# Contains data for initial conditions

test_options = ["test1-2", "noh", "test1-1"]

# diaphragm locations
dia_locs = {
    "test1-2": 0.5,
    "noh": 0.5,
    "test1-1": 0.3
}
# left and right states
rho_left = {
    "test1-2": 1.0,
    "noh": 1.0,
    "test1-1": 1.0
}
rho_right = {
    "test1-2": 1.0,
    "noh": 1.0,
    "test1-1": 0.125
}
p_left = {
    "test1-2": 0.4,
    "noh": 1e-6,
    "test1-1": 1.0
}
p_right = {
    "test1-2": 0.4,
    "noh": 1e-6,
    "test1-1": 0.1
}
u_left = {
    "test1-2": -2.0,
    "noh": 1.0,
    "test1-1": 0.75
}
u_right = {
    "test1-2": 2.0,
    "noh": -1.0,
    "test1-1": 0.0
}
# end times
end_times = {
    "test1-2": 0.15,
    "noh": 1,
    "test1-1": 0.2
}
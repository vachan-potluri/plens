# Contains data for initial conditions

test_options = [
    "noh",
    "test1-1",
    "test1-2",
    "test1-3"
]

# diaphragm locations
dia_locs = {
    "noh": 0.5,
    "test1-1": 0.3,
    "test1-2": 0.5,
    "test1-3": 0.5
}
# left and right states
rho_left = {
    "noh": 1.0,
    "test1-1": 1.0,
    "test1-2": 1.0,
    "test1-3": 1.0,
}
rho_right = {
    "noh": 1.0,
    "test1-1": 0.125,
    "test1-2": 1.0,
    "test1-3": 1.0
}
p_left = {
    "noh": 1e-6,
    "test1-1": 1.0,
    "test1-2": 0.4,
    "test1-3": 1000.0
}
p_right = {
    "noh": 1e-6,
    "test1-1": 0.1,
    "test1-2": 0.4,
    "test1-3": 0.01
}
u_left = {
    "noh": 1.0,
    "test1-1": 0.75,
    "test1-2": -2.0,
    "test1-3": 0.0
}
u_right = {
    "noh": -1.0,
    "test1-1": 0.0,
    "test1-2": 2.0,
    "test1-3": 0.0
}
# end times
end_times = {
    "noh": 1,
    "test1-1": 0.2,
    "test1-2": 0.15,
    "test1-3": 0.012
}
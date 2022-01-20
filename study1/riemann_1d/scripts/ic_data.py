# Contains data for initial conditions

test_options = [
    "noh",
    "test1-1",
    "test1-2",
    "test1-3",
    "test1-4"
]

# diaphragm locations
dia_locs = {
    "noh": 0.5,
    "test1-1": 0.3,
    "test1-2": 0.5,
    "test1-3": 0.5,
    "test1-4": 0.3 # 0.4 was given in Toro's book, but using 0.3
}
# left and right states
rho_left = {
    "noh": 1.0,
    "test1-1": 1.0,
    "test1-2": 1.0,
    "test1-3": 1.0,
    "test1-4": 5.99924
}
rho_right = {
    "noh": 1.0,
    "test1-1": 0.125,
    "test1-2": 1.0,
    "test1-3": 1.0,
    "test1-4": 5.99242
}
p_left = {
    "noh": 1e-6,
    "test1-1": 1.0,
    "test1-2": 0.4,
    "test1-3": 1000.0,
    "test1-4": 460.894
}
p_right = {
    "noh": 1e-6,
    "test1-1": 0.1,
    "test1-2": 0.4,
    "test1-3": 0.01,
    "test1-4": 46.095
}
u_left = {
    "noh": 1.0,
    "test1-1": 0.75,
    "test1-2": -2.0,
    "test1-3": 0.0,
    "test1-4": 19.5975
}
u_right = {
    "noh": -1.0,
    "test1-1": 0.0,
    "test1-2": 2.0,
    "test1-3": 0.0,
    "test1-4": -6.19633
}
# end times
end_times = {
    "noh": 1,
    "test1-1": 0.2,
    "test1-2": 0.15,
    "test1-3": 0.012,
    "test1-4": 0.035
}
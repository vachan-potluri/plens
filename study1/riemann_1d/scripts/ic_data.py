# Contains data for initial conditions

test_options = [
    "noh",
    "test1-1",
    "test1-2",
    "test1-3",
    "test1-4",
    "test1-5"
]

# diaphragm locations
dia_locs = {
    "noh": 0.5,
    "test1-1": 0.3,
    "test1-2": 0.5,
    "test1-3": 0.5,
    "test1-4": 0.3, # 0.4 was given in Toro's book, but using 0.3
    "test1-5": 0.8
}
# left and right states
rho_left = {
    "noh": 1.0,
    "test1-1": 1.0,
    "test1-2": 1.0,
    "test1-3": 1.0,
    "test1-4": 5.99924,
    "test1-5": 1.0
}
rho_right = {
    "noh": 1.0,
    "test1-1": 0.125,
    "test1-2": 1.0,
    "test1-3": 1.0,
    "test1-4": 5.99242,
    "test1-5": 1.0
}
p_left = {
    "noh": 1e-6,
    "test1-1": 1.0,
    "test1-2": 0.4,
    "test1-3": 1000.0,
    "test1-4": 460.894,
    "test1-5": 1000.0
}
p_right = {
    "noh": 1e-6,
    "test1-1": 0.1,
    "test1-2": 0.4,
    "test1-3": 0.01,
    "test1-4": 46.095,
    "test1-5": 0.01
}
u_left = {
    "noh": 1.0,
    "test1-1": 0.75,
    "test1-2": -2.0,
    "test1-3": 0.0,
    "test1-4": 19.5975,
    "test1-5": -19.59745
}
u_right = {
    "noh": -1.0,
    "test1-1": 0.0,
    "test1-2": 2.0,
    "test1-3": 0.0,
    "test1-4": -6.19633,
    "test1-5": -19.59745
}
# end times
end_times = {
    "noh": 1,
    "test1-1": 0.2,
    "test1-2": 0.15,
    "test1-3": 0.012,
    "test1-4": 0.035,
    "test1-5": 0.012
}
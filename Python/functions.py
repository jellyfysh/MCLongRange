import math

sigma = 1.0
epsilon = 1.0 / 0.46


# periodic distance between points a and b in a box of arbitrary size in any dimension
def per_dist(a, b, size):
    x_distances = [min(abs(a[j] - b[j]), size - abs(a[j] - b[j])) for j in range(len(a))]
    return math.sqrt(sum([x ** 2 for x in x_distances]))


# Lennard-Jones potential
def u_lj(r):
    return 4.0 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)


# derivative of u_lj
def du_lj(r):
    return 4 * epsilon * ((-12 * sigma ** 12) * r ** (-13) + 6 * (sigma ** 6) * r ** (-7))


# r-displacement producing a variation delta_u in u_lj, starting from r0
def r_disp(r0, delta_u):
    if 1.0 - math.sqrt(1.0 + (u_lj(r0) + delta_u) / epsilon) == 0.0:
        return [1.0, math.inf]
    r1 = (0.5 * (1.0 - math.sqrt(1.0 + (u_lj(r0) + delta_u) / epsilon))) ** (-1 / 6)
    r2 = (0.5 * (1.0 + math.sqrt(1.0 + (u_lj(r0) + delta_u) / epsilon))) ** (-1 / 6)
    return [r1, r2] if r1.imag == 0.0 else [-1.0, r2]


# x-displacement relative to r_final, starting from initial_pos
def x_disp(initial_pos, target_particle, r_final, size):
    delta_x_0 = per_dist([target_particle[0]], [initial_pos[0]], size)
    delta_y = per_dist([initial_pos[1]], [target_particle[1]], size)
    delta_x_final = math.sqrt(r_final ** 2 - delta_y ** 2) if delta_y <= r_final else 0.0
    return abs(delta_x_final - delta_x_0)


# x-displacement of active_particle producing a variation delta_u in u_lj with respect to target_particle
def displacement_advanced(active_particle, target_particle, delta_u, size):
    if delta_u == 0.0:
        return 0.0
    r_min = per_dist([target_particle[0], active_particle[1]], target_particle, size)
    r_max = per_dist([(target_particle[0] + size / 2) % size, active_particle[1]], target_particle, size)
    r0 = per_dist(active_particle, target_particle, size)
    r_well = 2 ** (1 / 6) * sigma
    # variation of u_lj over a period
    if r_min > r_well or r_max < r_well:
        delta_u_period = abs(u_lj(r_min) - u_lj(r_max))
    else:
        delta_u_period = abs(u_lj(r_min) - u_lj(r_well)) + abs(u_lj(r_max) - u_lj(r_well))
    n_periods = int(delta_u / delta_u_period)
    du = delta_u_period * n_periods
    u_left = delta_u - du
    if u_left <= 0.0:
        return size * n_periods
    r_move = per_dist([(active_particle[0] + 10 ** (-3)) % size, active_particle[1]], target_particle, size)
    moving_away = 1 if r_move > r0 else 0  # direction of motion of the active particle
    x_left1 = x_left2 = x_left3 = x_left4 = 0.0
    if r_min >= r_well:
        if moving_away:
            r = max(r_disp(r0, u_left))
            if r0 <= r <= r_max:
                x_left1 = x_disp(active_particle, target_particle, r, size)
            else:
                x_left1 = x_disp(active_particle, target_particle, r_max, size)
                x_left2 = size / 2
                u_left = u_left - (u_lj(r_max) - u_lj(r0))
                r = max(r_disp(r_min, u_left))
                x_left3 = x_disp([target_particle[0], active_particle[1]], target_particle, r, size)
        if not moving_away:
            x_left1 = x_disp(active_particle, target_particle, r_min, size)
            r = max(r_disp(r_min, u_left))
            x_left2 = x_disp([target_particle[0], active_particle[1]], target_particle, r, size)
    if r_max <= r_well:
        # this case never occurs if the box is big enough
        if moving_away:
            x_left1 = x_disp(active_particle, target_particle, r_max, size)
            r = min(r_disp(r_max, u_left)) if min(r_disp(r_max, u_left)) > 0.0 else max(r_disp(r_max, u_left))
            x_left2 = x_disp([(target_particle[0] + size / 2) % size, active_particle[1]], target_particle, r, size)
        if not moving_away:
            r = min(r_disp(r0, u_left)) if min(r_disp(r0, u_left)) > 0.0 else max(r_disp(r0, u_left))
            if r_min <= r <= r0:
                x_left1 = x_disp(active_particle, target_particle, r, size)
            else:
                x_left1 = x_disp(active_particle, target_particle, r_min, size)
                x_left2 = size / 2
                u_left = u_left - (u_lj(r_min) - u_lj(r0))
                r = min(r_disp(r_max, u_left)) if min(r_disp(r_max, u_left)) > 0 else max(r_disp(r_max, u_left))
                x_left3 = x_disp([(target_particle[0] + size / 2) % size, active_particle[1]], target_particle, r, size)
    if r_min < r_well < r_max:
        delta_y = per_dist([active_particle[1]], [target_particle[1]], size)
        x_well = [(target_particle[0] + math.sqrt(r_well ** 2 - delta_y ** 2)) % size, active_particle[1]]
        if r0 >= r_well:
            if moving_away:
                r = max(r_disp(r0, u_left))
                if r0 <= r <= r_max:
                    x_left1 = x_disp(active_particle, target_particle, r, size)
                else:
                    x_left1 = x_disp(active_particle, target_particle, r_max, size)
                    x_left2 = x_disp([(target_particle[0] + size / 2) % size, active_particle[1]],
                                     target_particle, r_well, size)
                    u_left = u_left - (u_lj(r_max) - u_lj(r0))
                    r = min(r_disp(r_well, u_left)) if min(r_disp(r_well, u_left)) > 0 else max(r_disp(r_well, u_left))
                    if r_min < r < r_well:
                        x_left3 = x_disp(x_well, target_particle, r, size)
                    else:
                        x_left3 = 2 * x_disp(x_well, target_particle, r_min, size)
                        u_left = u_left - (u_lj(r_min) - u_lj(r_well))
                        r = max(r_disp(r_well, u_left))
                        x_left4 = x_disp(x_well, target_particle, r, size)
            if not moving_away:
                x_left1 = x_disp(active_particle, target_particle, r_well, size)
                r = min(r_disp(r_well, u_left)) if min(r_disp(r_well, u_left)) > 0 else max(r_disp(r_well, u_left))
                if r_min <= r <= r_well:
                    x_left2 = x_disp(x_well, target_particle, r, size)
                else:
                    x_left2 = 2 * x_disp(x_well, target_particle, r_min, size)
                    u_left = u_left - (u_lj(r_min) - u_lj(r_well))
                    r = max(r_disp(r_well, u_left))
                    x_left3 = x_disp(x_well, target_particle, r, size)
        if r0 < r_well:
            if moving_away:
                x_left1 = x_disp(active_particle, target_particle, r_well, size)
                r = max(r_disp(r_well, u_left))
                if r_well <= r <= r_max:
                    x_left2 = x_disp(x_well, target_particle, r, size)
                else:
                    x_left2 = 2 * x_disp(x_well, target_particle, r_max, size)
                    u_left = u_left - (u_lj(r_max) - u_lj(r_well))
                    r = min(r_disp(r_well, u_left)) if min(r_disp(r_well, u_left)) > 0 else max(r_disp(r_well, u_left))
                    x_left3 = x_disp(x_well, target_particle, r, size)
            if not moving_away:
                r = min(r_disp(r0, u_left)) if min(r_disp(r0, u_left)) > 0 else max(r_disp(r0, u_left))
                if r_min <= r <= r0:
                    x_left1 = x_disp(active_particle, target_particle, r, size)
                else:
                    x_left1 = x_disp(active_particle, target_particle, r_min, size)
                    x_left2 = x_disp([target_particle[0], active_particle[1]], target_particle, r_well, size)
                    u_left = u_left - (u_lj(r_min) - u_lj(r0))
                    r = max(r_disp(r_well, u_left))
                    if r_well <= r <= r_max:
                        x_left3 = x_disp(x_well, target_particle, r, size)
                    else:
                        x_left3 = 2 * x_disp(x_well, target_particle, r_max, size)
                        u_left = u_left - (u_lj(r_max) - u_lj(r_well))
                        r = min(r_disp(r_well, u_left))
                        x_left4 = x_disp(x_well, target_particle, r, size)
    return size * n_periods + (x_left1 + x_left2 + x_left3 + x_left4)

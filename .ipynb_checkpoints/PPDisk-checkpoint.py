


class PPDisk():
    '''

    '''

    G = const.G.value
    M_s = const.M_sun.value
    AU = const.au.value

    def __init__(self, m_star, d_pc, sys_vel, inc, pa):
        '''
        Physical properties of protoplanetary disk
        '''
        self.m_star = m_star
        self.d_pc = d_pc
        self.sys_vel = sys_vel
        self.inc = inc
        self.pa = pa

    def

from ..imports import *

__all__ = ['tidy_radius_ticks']

def tidy_radius_ticks():
    '''
    Call plt.yticks(*tidy_radius_ticks()) to get nicer ticks for
    a log scale on a vertical radius axis (or likewise for horizontal).
    '''
    t = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30]
    s = ['1', '2', '3', '4', '5', ' ', ' ', ' ', ' ', '10', '20', '30']

    return t, s

def fusswithticks(panel):
    '''
    FIXME -- can probably remove?
    '''
    plt.sca(panel.ax)
    t = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30]
    s = ['1', '2', '3', '4', ' ', ' ', ' ', ' ', ' ', '10', '20']
    plt.yticks(t, s)


    # mass on the
    t = [.7, .8, .9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 0.5]
    s = ['{0}'.format(v) for v in t]
    s[0] = ' '
    s[1] = ' '
    s[2] = ' '
    s[8] = ' '
    s[9] = ' '
    s[10] = ' '
    s[11] = ' '


    plt.xticks(t, s)

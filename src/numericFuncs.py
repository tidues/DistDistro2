
# numeric gate function
def eta(a, b, x):
    if x >= a and x <=b:
        return 1
    else:
        return 0

# numeric threshold function
def theta(a, b, x):
    if x < a:
        return a
    elif x >b:
        return b
    else:
        return x


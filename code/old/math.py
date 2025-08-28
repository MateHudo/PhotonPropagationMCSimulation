def sum(x,y):
    print("add numbers {} and {}".format(x,y))
    return x+y
print(sum(1,2))

import matplotlib.pyplot as plt

## new code block - define subtraction
def subtract(x,y):
    print('subtract {} with {}'.format(x,y))
    return x-y
# example
print(subtract(3,2))


# Sample data
x = [1, 2, 3, 4, 5]
y = [1, 4, 9, 16, 25]

# Create a plot
plt.plot(x, y)
# Add titles and labels
plt.title('Simple Plot')
plt.xlabel('time')
plt.ylabel('loaction')
# Show the plot
plt.show()


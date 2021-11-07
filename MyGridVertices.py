import pygame, sys
from pygame.locals import *

# Point and Line class definitions
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Line:
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2

    def getStartPoint(self):
        return (self.p1.x, self.p1.y)

    def getEndPoint(self):
        return (self.p2.x, self.p2.y)

# Verifies if 3 points are colinear
def onLine(p, q, r):
    if ( (q.x <= max(p.x, r.x)) and (q.x >= min(p.x, r.x)) and
           (q.y <= max(p.y, r.y)) and (q.y >= min(p.y, r.y))):
        return True
    return False

# Finds orientation of 3 points, clockwise, counterclockwise or collinear
def orientation(p, q, r):
    # to find the orientation of 3 points (p, q, r)
    # function returns the following values:
    # 0 : Collinear points
    # 1 : Clockwise points
    # 2 : Counterclockwise

    val = (float(q.y - p.y) * (r.x - q.x)) - (float(q.x - p.x) * (r.y - q.y))
    if (val > 0):

        # Clockwise orientation
        return 1
    elif (val < 0):

        # Counterclockwise orientation
        return 2
    else:

        # Collinear orientation
        return 0


# The main function that returns true if
# the line segment 'p1q1' and 'p2q2' intersect.
def doIntersect(line1, line2):
    # Find the 4 orientations required for
    # the general and special cases

    p1 = line1.p1
    q1 = line1.p2
    p2 = line2.p1
    q2 = line2.p2

    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    # General case
    if ((o1 != o2) and (o3 != o4)):
        return True

    # Special Cases

    # p1 , q1 and p2 are collinear and p2 lies on segment p1q1
    if ((o1 == 0) and onLine(p1, p2, q1)):
        return True

    # p1 , q1 and q2 are collinear and q2 lies on segment p1q1
    if ((o2 == 0) and onLine(p1, q2, q1)):
        return True

    # p2 , q2 and p1 are collinear and p1 lies on segment p2q2
    if ((o3 == 0) and onLine(p2, p1, q2)):
        return True

    # p2 , q2 and q1 are collinear and q1 lies on segment p2q2
    if ((o4 == 0) and onLine(p2, q1, q2)):
        return True

    # If none of the cases
    return False

# Simple method to draw a straight line after validation
def drawLine(surface, newLine):
    start = newLine.getStartPoint()
    end = newLine.getEndPoint()
    pygame.draw.line(surface, black, start, end)
    pygame.display.update()

# Return true if line does not intersect, false if it does
def validateLine(newLine, lineList):
    intersect = False

    if(len(lineList) == 0):
        return True
    else:
        for line in lineList:
            if (line == lineList[len(lineList) - 1]):
                return True
            intersect = doIntersect(line, newLine)
            if(intersect):
                return False

    return True

# Validation, called only when trying to exit program
def validateLastLine(newLine, lineList):
    if (len(lineList) == 0):
        return True
    else:
        tempList = [line for line in lineList]
        del tempList[0]
        for line in tempList:
            if (line == tempList[len(tempList) - 1]):
                return True
            intersect = doIntersect(line, newLine)
            if (intersect):
                return False

    return True

# Method that controls program flow
def continueShape(x, y):
    newPoint = Point(x, y)
    if (len(points) == 0):
        points.append(newPoint)
    elif (len(points) == 1):
        newLine = Line(points[len(points) - 1], newPoint)
        points.append(newPoint)
        lines.append(newLine)
        drawLine(grid, newLine)
    else:
        newLine = Line(points[len(points) - 1], newPoint)
        if (validateLine(newLine, lines)):
            points.append(newPoint)
            lines.append(newLine)
            drawLine(grid, newLine)

# Called at the end to finish the shape when user exits
def closeShape(x, y):
    newPoint = Point(x, y)
    if (len(points) == 0):
        points.append(newPoint)
    elif (len(points) == 1):
        newLine = Line(points[len(points) - 1], newPoint)
        points.append(newPoint)
        lines.append(newLine)
        drawLine(grid, newLine)
    else:
        newLine = Line(points[len(points) - 1], newPoint)
        if (validateLastLine(newLine, lines)):
            points.append(newPoint)
            lines.append(newLine)
            drawLine(grid, newLine)

# Returns correctly formatted output: List of (x,y) coordinates of all points
def output(pointList):
    pointListReversed = pointList[::-1]
    output = [[float(point.x), float(point.y)] for point in pointListReversed]
    return output

# Color definitions
black = 0, 0, 0
white = 255, 255, 255
red = 255, 0, 0

# Define grid size
gridWidth, gridHeight = 1000, 800
gridSize = gridWidth, gridHeight

# Setting the grid
grid = pygame.display.set_mode(gridSize)

# Getting the Clock object
clock = pygame.time.Clock()

# Settomg a title to the grid
pygame.display.set_caption("Grid")

# Defining lists for all points and all lines
points = []
lines = []

def programLoop():
    fps_limit = 60
    run_me = True
    while run_me:
        # Limit the framerate
        clock.tick(fps_limit)

        # Clear the grid by filling all white
        if (len(points) == 0):
            grid.fill(white)

        # Event handler
        for event in pygame.event.get():
            if event.type == MOUSEBUTTONDOWN:
                x, y = pygame.mouse.get_pos()
                continueShape(x, y)
            if event.type == pygame.QUIT:
                numOfPoints = len(points)
                closeShape(points[0].x, points[0].y)
                if (numOfPoints == len(points)):
                    continue

                # Final output, must decide what to do with it
                output(points)
                run_me = False

        # Display everything
        pygame.display.flip()

programLoop()


# Quit the display
pygame.quit()
sys.exit()
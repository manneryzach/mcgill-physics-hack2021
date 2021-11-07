import pygame, sys
from pygame.locals import *
import math

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

class Grid:

    def __init__(self, grid, gridWidth, gridHeight, squareWidth):
        self.gridWidth = gridWidth
        self.gridHeight = gridHeight
        self.squareWidth = squareWidth
        self.numCols = int(gridWidth / squareWidth)
        self.numRows = int(gridHeight / squareWidth)
        self.grid = grid
        self.gridArr = [[0 for i in range(self.numCols)] for j in range(self.numRows)]

    def select(self, x, y):
        xIndex = math.floor(x / squareWidth)
        yIndex = math.floor(y / squareWidth)
        print(type(self.gridArr[yIndex][xIndex]))
        self.gridArr[yIndex][xIndex].changeColor(red)

    def update(self):
        for currentY in range(0, self.gridHeight, self.squareWidth):
            for currentX in range(0, gridWidth, self.squareWidth):
                self.gridArr[int(currentY / self.squareWidth)][int(currentX / self.squareWidth)].display()

    def display(self):
        # Fill grid with square objects. Always divide by squareWidth for indices of list
        for currentY in range(0, self.gridHeight, self.squareWidth):
            for currentX in range(0, gridWidth, self.squareWidth):
                self.gridArr[int(currentY / self.squareWidth)][int(currentX / self.squareWidth)] = GridSquare(currentX, currentY,
                                                                                              self.squareWidth, white, 1)
                self.gridArr[int(currentY / self.squareWidth)][int(currentX / self.squareWidth)].display()

        # Draw horizontal grid lines.
        for currentY in range(0, self.gridHeight, self.squareWidth):
            pygame.draw.line(grid, black, (0, currentY), (self.gridWidth, currentY))
            if (currentY == self.gridHeight - self.squareWidth):
                pygame.draw.line(grid, black, (0, currentY + self.squareWidth), (self.gridWidth, currentY + self.squareWidth))

        # Draw vertical grid lines.
        for currentX in range(0, self.gridWidth, self.squareWidth):
            pygame.draw.line(grid, black, (currentX, 0), (currentX, self.gridHeight))
            if (currentX == self.gridWidth - self.squareWidth):
                pygame.draw.line(grid, black, (currentX + self.squareWidth, 0), (currentX + self.squareWidth, self.gridHeight))

class GridSquare:

    def __init__(self, x, y, size, color=(255,255,255), lineWidth=1):
        self.x = x
        self.y = y
        self.size = size
        self.color = color
        self.lineWidth = lineWidth

    def display(self):
        pygame.draw.polygon(grid, self.color, [(self.x, self.y), (self.x + self.size, self.y), (self.x + self.size, self.y + self.size), (self.x, self.y + self.size)], width=0)

    def changeColor(self, color):
        self.color = color



testSquare = GridSquare(0, 20, 10, black, 1)

# Defining variables for fps and continued running
fps_limit = 60
run_me = True
while run_me:
    # Limit the framerate
    clock.tick(fps_limit)

    # Clear the screen by filling all white
    grid.fill(white)

    # Display Grid
    squareWidth = 20
    myGrid = Grid(grid, gridWidth, gridHeight, squareWidth)
    myGrid.display()

    for event in pygame.event.get():
        if event.type == MOUSEBUTTONDOWN:
            x, y = pygame.mouse.get_pos()
            print(x)
            print(y)
            myGrid.select(x, y)
            myGrid.update()
        if event.type == pygame.QUIT:
            run_me = False


    # Display everything
    pygame.display.flip()


# Quit the display
pygame.quit()
sys.exit()
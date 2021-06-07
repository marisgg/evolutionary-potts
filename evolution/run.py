import evolution

F = "src/runForagingModel.js"
A = "src/runActModel.js"

# This main is needed on Windows to prevent recursive subprocesses opening
if __name__ in '__main__':
    evolution.evolve(F, 1)
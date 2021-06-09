import evolution

F = "src/runForagingModel.js"
A = "src/runActModel.js"

# This main is needed on Windows to prevent recursive subprocesses opening
if __name__ in '__main__':
    seed = None
    if seed is None:
        print("Enter a starting seed:")
        seed = int(input())
    evolution.evolve(F, 30, seed)
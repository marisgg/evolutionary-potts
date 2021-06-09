import evolution
import os
F = "src/runForagingModel.js"
F1 = "src/runForagingModel-mediumchemotaxis.js"
F2 = "src/runForagingModel-nochemotaxis.js"
A = "src/runActModel.js"

# This main is needed on Windows to prevent recursive subprocesses opening
if __name__ in '__main__':
    print(f"Enter chemotaxis level (none, medium, high):")
    level = input()
    if level == "high":
        model = F
    elif level == "medium":
        model = F1
    elif level == "none":
        model = F2
    else:
        raise ValueError("Incorrect chemotaxis level string")
    if not os.path.exists(model):
        raise FileNotFoundError(model)
    seed = None
    if seed is None:
        print("Enter a starting seed:")
        seed = int(input())
    evolution.evolve(model, 30, seed)
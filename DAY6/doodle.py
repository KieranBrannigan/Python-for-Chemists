a = [1,2,3,4,5,6,7]
b = [1,2,3,4,5,6,7]

for val in a:
    for i in b:
        print(f"val={val}, i={i}")
        if i == val:
            print("inside the if statement , calling break")
            break
        else:
            print("else statement")
            

        print("outside the if else loop, calling break")
        break

    
    print("outside the for i in b loop")
        
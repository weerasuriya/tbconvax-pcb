## New 3 age-class demographic model

- In the conventional single age model
    - When calculating y = 2, step 1
        - The receiving grid is all zero
        - ages 0:max-1, i-1 contents go into receiving grid age 1:max, step i
        - all compartments of age=0, step i are empty and zerod
        - births are added to age=0, step i, but no ops happen to this group
    - y = 2, step = 2
        - receiving grid all zero
        - ages 0:max go to ages 0:max, step i-1 to step i
        - age 0 undergoes usual dynamic movement

- In the 3 age-class model
    - When calculating y = 2, step 1
        - The receiving grid is all zero
        -

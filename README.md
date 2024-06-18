# PQTD_MCTS
## Background
- Implementation of task dependency assignment.
- Prioritized task assignment approach based on MCTS algorithm.
## Structure and important parameters
- **test.cpp**: Contains the main function for testing.
- **daggen_commons.cpp**: Sets up task information and server information.
- **mcts.cpp**: Implements the PQTD+MCTS algorithm.

## Usage
- For Windows, you need to have Visual Studio 2019 installed.
- For Linux, you need to have the GCC compiler installed.
### Using in VS2019
1. Open Visual Studio 2019.

2. Clone the repository and open the solution file 'test' in VS2019.

3. In the Solution Explorer, right-click on the project and select "Build".

5. After building, locate the generated executable file and run `test.cpp`.

### Using in Linux
1. Clone the repository:
    ```bash
    git clone https://github.com/maleficent1997/PQTD_MCTS.git
    cd [the file path]/test_linux
    ```
2. Compile the code using g++:
    ```bash
    g++ test.cpp mcts.cpp daggen_commons.cpp -o test
    ```
3. Run the generated executable:
    ```bash
    ./test
    ```
## Contact
Email : caiqing19s@ict.ac.cn

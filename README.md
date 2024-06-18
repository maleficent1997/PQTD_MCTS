# PQTD_MCTS
## Background
The open source code of the paper 'Prioritized Assignment with Task Dependency in Collaborative Mobile Edge Computing'.
## Structure and important parameters

## Usage
- For Windows, you need to have Visual Studio 2019 installed.
- For Linux, you need to have the GCC compiler installed.
### Using in VS2019
1. Open Visual Studio 2019.

2. Clone the repository and open the solution file 'test' in in VS2019.

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

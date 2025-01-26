# Costas Array Generator

This repository provides tools to generate and visualize **Costas arrays**, mathematical structures with applications in radar, sonar, and communications.

## Features

- **Welch Construction** for fast generation when applicable (\(N+1\) is prime).
- **Backtracking Search** as a fallback for general \(N\).
- **Visualization** of arrays using `matplotlib`.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/splch/costas-array-generation.git
   cd costas-array-generation
   ```
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

```python
from costas import generate_costas_array, visualize_costas_array

N = 8
result = generate_costas_array(N)

print("Costas permutation:", result)
visualize_costas_array(result)
```

![image](https://github.com/user-attachments/assets/37c24180-7f08-4ed2-a8b5-0bd0f66284d4)

Run examples in the provided `main.ipynb`.

## License

This project is open-sourced under the [MIT License](https://opensource.org/licenses/MIT).

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description='Plot L2 error vs h from a CSV file.')
    parser.add_argument('filename', nargs='?', default='L2file.csv', help='Input CSV file (default: L2file.csv)')
    args = parser.parse_args()

    data = pd.read_csv(args.filename)

    h = np.array(data.h)
    L2 = np.array(data.L2_error)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Linear scale plot
    ax1.plot(h, L2, 'o-')
    ax1.set_title("L2-error versus h")
    ax1.set_xlabel("h")
    ax1.set_ylabel("L2 error")
    ax1.grid(True)

    # Log-log scale plot
    ax2.loglog(h, L2, 'o-')
    ax2.set_title("L2-error versus h (log-log)")
    ax2.set_xlabel("h")
    ax2.set_ylabel("L2 error")
    ax2.grid(True, which="both", ls="-")

    plt.tight_layout()
    plt.savefig("L2graph.png")
    print("Graph saved to L2graph.png")

if __name__ == "__main__":
    main()




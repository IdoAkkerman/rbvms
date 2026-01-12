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
    H1 = np.array(data.H1_error)

    fig,ax = plt.subplots(2, 2, figsize=(12, 10))

    # Linear scale plot
    ax[0,0].plot(h, L2, 'o-')
    ax[0,0].set_title("L2-error versus h")
    ax[0,0].set_xlabel("h")
    ax[0,0].set_ylabel("L2 error")
    ax[0,0].grid(True)

    # Log-log scale plot
    ax[0,1].loglog(h, L2, 'o-')
    ax[0,1].set_title("L2-error versus h (log-log)")
    ax[0,1].set_xlabel("h")
    ax[0,1].set_ylabel("L2 error")
    ax[0,1].grid(True, which="both", ls="-")

    # Linear scale plot
    ax[1,0].plot(h, H1, 'o-')
    ax[1,0].set_title("H1-error versus h")
    ax[1,0].set_xlabel("h")
    ax[1,0].set_ylabel("H1 error")
    ax[1,0].grid(True)

    # Log-log scale plot
    ax[1,1].loglog(h, H1, 'o-')
    ax[1,1].set_title("H1-error versus h (log-log)")
    ax[1,1].set_xlabel("h")
    ax[1,1].set_ylabel("H1 error")
    ax[1,1].grid(True, which="both", ls="-")

    plt.tight_layout()
    plt.savefig("L2graph.png")
    print("Graph saved to L2graph.png")

if __name__ == "__main__":
    main()




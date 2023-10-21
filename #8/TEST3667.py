import tkinter as tk
from tkinter import filedialog
import cv2
import numpy as np
import os


# Function to compute EHD descriptor
def compute_ehd_descriptor(image, num_bins=36):
    # Convert image to grayscale
    gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    
    # Compute gradients in the x and y directions
    grad_x = cv2.Sobel(gray_image, cv2.CV_64F, 1, 0, ksize=3)
    grad_y = cv2.Sobel(gray_image, cv2.CV_64F, 0, 1, ksize=3)
    
    # Compute magnitude and orientation of gradients
    magnitude = np.sqrt(grad_x**2 + grad_y**2)
    orientation = np.arctan2(grad_y, grad_x) * (180 / np.pi)
    
    # Ensure orientation values are positive
    orientation[orientation < 0] += 360
    
    # Create histogram bins
    bin_width = 360 / num_bins
    bins = np.arange(0, 360, bin_width)
    
    # Compute histogram of edge orientations
    hist, _ = np.histogram(orientation, bins=bins, range=(0, 360))
    
    # Normalize the histogram
    hist = hist.astype(float) / np.sum(hist)
    
    return hist

# Function to compute similarity/distance between descriptors
def compute_distance(descriptor1, descriptor2):
    return np.linalg.norm(descriptor1 - descriptor2)

# Function to search for similar images
def search_similar_images(query_image_path, folder_path, num_results=5):
    # Perform image search using the stored paths
    if not query_image_path or not folder_path:
        print("Please select the query image and folder.")
        return

    # Load query image
    print(query_image_path)
    query_image = cv2.imread(query_image_path)
    
    # Compute EHD descriptor for query image
    query_descriptor = compute_ehd_descriptor(query_image)
    
    # Initialize list to store image paths and distances
    image_paths = []
    distances = []
    
    # Iterate over images in the folder
    for filename in os.listdir(folder_path):
        image_path = os.path.join(folder_path, filename)
        # Skip non-image files
        if not cv2.haveImageReader(image_path):
            continue
        
        # Load and compute descriptor for each image in the folder
        image = cv2.imread(image_path)
        image_descriptor = compute_ehd_descriptor(image)
        
        # Compute distance between query descriptor and image descriptor
        distance = compute_distance(query_descriptor, image_descriptor)
        
        # Store image path and distance
        image_paths.append(image_path)
        distances.append(distance)
    
    # Sort image paths based on distances
    sorted_paths = [x for _, x in sorted(zip(distances, image_paths))]
    sorted_distances = sorted(distances)
    
    # Display the most similar images
    for i in range(num_results):
        image = cv2.imread(sorted_paths[i])
        cv2.imshow(f"Similar Image {i+1}", image)
    
    cv2.waitKey(0)
    cv2.destroyAllWindows()

# Function to open file dialog and select image file
def select_query_image():
    query_image_path.set(filedialog.askopenfilename(title="Select Query Image"))

# Function to open file dialog and select folder
def select_folder():
    folder_path.set(filedialog.askdirectory(title="Select Folder"))
    

# Function to handle browse button for query image
def browse_query_image():
    filename = filedialog.askopenfilename(initialdir="/", title="Select Query Image",
                                          filetypes=(("JPEG files", "*.jpg"), ("PNG files", "*.png")))
    query_entry.delete(0, tk.END)
    query_entry.insert(0, filename)


# Function to handle browse button for folder
def browse_folder():
    foldername = filedialog.askdirectory(initialdir="/", title="Select Image Folder")
    folder_entry.delete(0, tk.END)
    folder_entry.insert(0, foldername)

# Function to handle run button
def run_search():
    num_results = 5
    query_image_path1 = query_image_path.get()
    folder1 = folder_path.get()
    print(query_image_path.get())
    print(folder_path.get())
    
    # Validate the inputs
    if not query_image_path or not folder_path or not num_results:
        print("Please select the query image, folder, and enter the number of results.")
        return

    try:
        num_results = int(num_results)
    except ValueError:
        print("Please enter a valid number of results.")
        return

    search_similar_images(query_image_path1, folder1, num_results)

# Create the main window
window = tk.Tk()
window.title("Image Search")
window.geometry("400x250")

# Create StringVar objects
query_image_path = tk.StringVar()
folder_path = tk.StringVar()
number_results= tk.StringVar()

# Create labels
query_label = tk.Label(window, text="Query Image:")
query_label.pack()
folder_label = tk.Label(window, text="Folder:")
folder_label.pack()
results_label = tk.Label(window, text="Number of Results:")
results_label.pack()

# Create entry fields
query_entry = tk.Entry(window)
query_entry.pack()
folder_entry = tk.Entry(window)
folder_entry.pack()
results_entry = tk.Entry(window)
results_entry.pack()

# Create buttons
query_button = tk.Button(window, text="Select Query Image", command=select_query_image)
print(query_image_path)
query_button.pack()
folder_button = tk.Button(window, text="Select Folder", command=select_folder)
folder_button.pack()
results_entry = tk.Entry(window, textvariable=number_results)
results_entry.pack()
run_button = tk.Button(window, text="Run", command=run_search)
run_button.pack()

# Create error label
error_label = tk.Label(window, fg="red")
error_label.pack()

# Start the main loop
window.mainloop()

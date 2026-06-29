import cv2 as cv
import numpy as np
import os
import tensorflow as tf
from keras.models import load_model
from keras.applications.resnet_v2 import preprocess_input

# Change as needed
MODEL_PATH = 'selectedModel.h5'
IMG_SIZE = 224
CLASS_NAMES = [
    "BacterialBlights",
    "Healthy",
    "Mosaic",
    "RedRot",
    "Rust",
    "Yellow"
  ] 

def square_resize_image(img_input, target_size):
    # Standarization
    img_out = img_input
    h, w = img_out.shape[:2]
    
    # Square-shaping the image
    top, bottom, left, right = 0, 0, 0, 0
    if w >= h:
        top = (w - h) // 2
        bottom = (w - h) - top
    else:
        left = (h - w) // 2
        right = (h - w) - left
        
    if top > 0 or bottom > 0 or left > 0 or right > 0:
        img_out = cv.copyMakeBorder(img_out, top, bottom, left, right, cv.BORDER_CONSTANT, value=(0, 0, 0))
    
    # Resize
    img_out = cv.resize(img_out, (target_size, target_size), interpolation=cv.INTER_AREA)
    return img_out

def load_and_preprocess_image(image_path):
    img = cv.imread(image_path)
    if img is None:
        print(f"The path {image_path} couldn't be read.")
        return None

    # Sizing
    img = square_resize_image(img, IMG_SIZE)

    # BGR -> RGB
    img = cv.cvtColor(img, cv.COLOR_BGR2RGB)
    img = img.astype('float32')

    # Adding an extra dimension (1, 224, 224, 3)
    img_expanded = np.expand_dims(img, axis=0)
    
    return img_expanded

def main():
    if not os.path.exists(MODEL_PATH):
        print(f"There's no model at {MODEL_PATH}")
        return
    
    try:
        model = load_model(MODEL_PATH, custom_objects={'preprocess_input': preprocess_input})
    except Exception as e:
        print(f"\Error while loading the model:\n{e}")
        return

    while True:
        image_path = input("\Introduce the image's route or type 'q' to cancel: ").strip()
        image_path = image_path.replace('"', '').replace("'", "")

        if image_path.lower() == 'q':
            break

        if not os.path.exists(image_path):
            print("This route doesn't exist.")
            continue

        processed_img = load_and_preprocess_image(image_path)
        
        if processed_img is not None:
            # Predicting
            predictions = model.predict(processed_img, verbose=0)
            
            # Translating
            predicted_class_idx = np.argmax(predictions[0])
            confidence = np.max(predictions[0]) * 100

            if predicted_class_idx < len(CLASS_NAMES):
                class_name = CLASS_NAMES[predicted_class_idx]
            else:
                class_name = f"Class {predicted_class_idx} (Unknown)"

            print(f"------------------------------------------------")
            print(f"Prediction: {class_name}")
            print(f"Confidence:  {confidence:.2f}%")
            print(f"------------------------------------------------")

if __name__ == "__main__":
    main()
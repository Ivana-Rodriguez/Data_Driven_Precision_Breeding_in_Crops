from settings import *

import cv2 as cv

import os

import random



# --- CONFIGURACIÓN DE SPLIT (80% Train, 10% Val, 10% Test) ---
TRAIN_RATIO = 0.8
VAL_RATIO = 0.1

# -----------------------------------------------------------
def square_resize_image(img_input):
    img_out = img_input
    if not stretch:
        w, h = img_out.shape[1], img_out.shape[0]
        top, bottom, left, right = 0, 0, 0, 0
        if w >= h:
            top = (w - h) // 2
            bottom = (w - h) - top
        else:
            left = (h - w) // 2
            right = (h - w) - left
        img_out = cv.copyMakeBorder(img_out, top, bottom, left, right, cv.BORDER_CONSTANT, value=(0, 0, 0))
    img_out = cv.resize(img_out, (image_size, image_size), interpolation=cv.INTER_AREA)
    return img_out

print(f"--- INICIANDO PROCESO ---")
print(f"Buscando clases en: {original_data_folder}")

if not os.path.exists(original_data_folder):
    print(f"ERROR: No existe la ruta {original_data_folder}")
    exit()

classes = []
for d in os.listdir(original_data_folder):
    path_completo = os.path.join(original_data_folder, d)

    # Verificamos que sea carpeta y que no sea la carpeta 'testing' si no la quieres
    if os.path.isdir(path_completo) and d != 'testing' and not d.startswith('.'):
        classes.append(d)
print(f"Clases encontradas ({len(classes)}): {classes}")
if len(classes) == 0:
    print("ERROR: Sigue sin encontrar carpetas. Revisa que 'root_folder' en settings.py sea correcta.")
    exit()

for cls in classes:
    cls_root_folder = os.path.join(original_data_folder, cls)

    # 2. RECOLECCIÓN RECURSIVA (Busca en todas las subcarpetas)
    all_image_paths = []
    for root, dirs, files in os.walk(cls_root_folder):
        for file in files:
            if file.lower().endswith(('.jpg', '.jpeg', '.png', '.tif', '.tiff', '.bmp')):
                full_path = os.path.join(root, file)
                all_image_paths.append(full_path)

    if not all_image_paths:
        print(f"AVISO: La clase '{cls}' existe pero está VACÍA.")
        continue
        
    # 3. MEZCLAR Y DIVIDIR
    random.shuffle(all_image_paths)
    total = len(all_image_paths)
    train_end = int(total * TRAIN_RATIO)
    val_end = int(total * (TRAIN_RATIO + VAL_RATIO))
    train_list = all_image_paths[:train_end]
    val_list = all_image_paths[train_end:val_end]
    test_list = all_image_paths[val_end:]
    print(f"Procesando '{cls}': Total {total} -> Train:{len(train_list)} | Val:{len(val_list)} | Test:{len(test_list)}")

    # 4. CREAR CARPETAS DE DESTINO
    train_cls_dir = os.path.join(train_folder, cls)
    os.makedirs(train_cls_dir, exist_ok=True)
    val_cls_dir = os.path.join(val_folder, cls)
    os.makedirs(val_cls_dir, exist_ok=True)
    test_cls_dir = os.path.join(test_folder, cls)
    os.makedirs(test_cls_dir, exist_ok=True)

    # FUNCIÓN DE GUARDADO
    def save_subset(path_list, target_dir):
        for full_path in path_list:
            img = cv.imread(full_path)
            if img is not None:
                img = square_resize_image(img)

                # Nombre base + extensión jpg
                filename = os.path.basename(full_path)
                name_no_ext = os.path.splitext(filename)[0]
                save_path = os.path.join(target_dir, name_no_ext + '.jpg')

                # Evitar sobrescribir si hay nombres repetidos
                counter = 1
                while os.path.exists(save_path):
                    save_path = os.path.join(target_dir, f"{name_no_ext}_{counter}.jpg")
                    counter += 1

                cv.imwrite(save_path, img)

    save_subset(train_list, train_cls_dir)
    save_subset(val_list, val_cls_dir)
    save_subset(test_list, test_cls_dir)
print("\n--- ¡LISTO! Revisa tus carpetas de data/train_... ---")
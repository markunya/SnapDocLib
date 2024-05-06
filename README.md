# SnapDocLib

SnapDocLib - это C++ библиотека для обработки изображений с использованием OpenCV.

## Возможности

- Обнаружение углов (DetectCorners)
- Нормализация изображений (Normalize)

## Сборка и установка

### Шаг 1: Установите OpenCV
Для установки OpenCV выполните следующие шаги:

1. Перейдите на [официальный сайт OpenCV](https://opencv.org/releases/) и скачайте версию 4.0 или выше.
2. Следуйте [инструкции по установке](https://docs.opencv.org/master/df/d65/tutorial_table_of_content_introduction.html) для вашей операционной системы (Windows, Linux или macOS).

### Шаг 2: Клонируйте репозиторий
Клонируйте репозиторий на локальную машину:
```bash
git clone https://github.com/markunya/SnapDocLib.git
cd SnapDocLib
```

### Шаг 3: Изменения в CMakeLists.txt
Измените путь в set(OpenCV_DIR "/home/xubuntu/OpenCV-android-sdk/sdk/native/jni"), на свой локальный. Без этой строчки у меня были проблемы со сборкой

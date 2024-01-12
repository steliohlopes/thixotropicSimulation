from PIL import Image

def overlay_images_with_white_background(image_path1, image_path2, output_path):
    # Abre as imagens
    img1 = Image.open(image_path1).convert("RGBA")
    img2 = Image.open(image_path2).convert("RGBA")

    # Remove o fundo branco
    img1 = remove_white_background(img1)
    img2 = remove_white_background(img2)

    # Cria uma nova imagem branca do mesmo tamanho que as originais
    white_background = Image.new("RGBA", img1.size, (255, 255, 255, 255))

    # Sobreposição das imagens sobre o fundo branco
    result = Image.alpha_composite(white_background, img1)
    result = Image.alpha_composite(result, img2)

    # Salva a imagem resultante
    result.save(output_path, "PNG")

def remove_white_background(image):
    # Obtém os dados dos pixels
    data = image.getdata()

    # Converte os pixels brancos (255, 255, 255, 255) para transparentes (0, 0, 0, 0)
    new_data = [(r, g, b, 0) if (r, g, b) == (255, 255, 255) else (r, g, b, a) for (r, g, b, a) in data]

    # Atualiza a imagem com os novos dados
    image.putdata(new_data)

    return image

# Exemplo de uso
overlay_images_with_white_background("PreProcessing/PipeFlow2D/SMDResult.png", "PreProcessing/PipeFlow2D/ThixotropicResultTa0.0001.png", "PreProcessing/PipeFlow2D/ThixotropicResultTa0.0001.png")

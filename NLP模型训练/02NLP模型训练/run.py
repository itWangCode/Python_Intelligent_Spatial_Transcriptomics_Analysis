from biosentvec_toolkit import BioSentVecToolkit

# 初始化（自动检测GPU/CPU）
toolkit = BioSentVecToolkit()

# 训练模型
toolkit.train_word_embeddings(
    corpus_file='medical_texts.txt',
    vector_size=200,
    epochs=10
)

# 计算相似度
sim = toolkit.compute_similarity(
    "Patient has diabetes",
    "Diabetes diagnosed"
)
print(f"相似度: {sim:.4f}")
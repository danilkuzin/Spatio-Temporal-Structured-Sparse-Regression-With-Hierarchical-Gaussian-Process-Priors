function result = NMSE(est_matrix, true_matrix)

result = trace((est_matrix - true_matrix)' * (est_matrix - true_matrix)) / trace(true_matrix' * true_matrix);

end
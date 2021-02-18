from unittest import TestCase

import os
import numpy as np
from data.test_data.test_parsers.config import test_parsers_folder_path
from pythhon.mesh.parsers import MeshFileContainer


class TestMeshFileContainer(TestCase):
    def test_vertices(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        p = MeshFileContainer(geof_file_path)
        vertices = np.array(
            [
                [0.0, 0.0],
                [0.5, 0.0],
                [1.0, 0.0],
                [0.0, 0.5],
                [0.5, 0.5],
                [1.0, 0.5],
                [0.0, 1.0],
                [0.5, 1.0],
                [1.0, 1.0],
            ]
        )
        bool_res_not = not (p.vertices - vertices).all()
        self.assertTrue(bool_res_not)

    def test_vertices_weights_cell_wise(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        p = MeshFileContainer(geof_file_path)
        vertices_weights_cell_wise_check = np.array(
            [1.0, 3.0, 2.0, 3.0, 6.0, 3.0, 2.0, 3.0, 1.0]
        )
        bool_res_not = not (
            p.vertices_weights_cell_wise - vertices_weights_cell_wise_check
        ).all()
        self.assertTrue(bool_res_not)

    def test_vertices_weights_face_wise(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        p = MeshFileContainer(geof_file_path)
        vertices_weights_face_wise = np.array(
            [2.0, 4.0, 3.0, 4.0, 6.0, 4.0, 3.0, 4.0, 2.0]
        )
        bool_res_not = not (
            p.vertices_weights_face_wise - vertices_weights_face_wise
        ).all()
        self.assertTrue(bool_res_not)

    def test_cells_vertices_connectivity_matrix(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        p = MeshFileContainer(geof_file_path)
        cells_vertices_connectivity_matrix = [
            np.array([0, 1, 3], dtype=np.uint8),
            np.array([1, 4, 3], dtype=np.uint8),
            np.array([3, 4, 6], dtype=np.uint8),
            np.array([4, 7, 6], dtype=np.uint8),
            np.array([1, 2, 4], dtype=np.uint8),
            np.array([2, 5, 4], dtype=np.uint8),
            np.array([4, 5, 7], dtype=np.uint8),
            np.array([5, 8, 7], dtype=np.uint8),
        ]
        bool_res_not = True
        for i in range(len(cells_vertices_connectivity_matrix)):
            if (
                p.cells_vertices_connectivity_matrix[i]
                - cells_vertices_connectivity_matrix[i]
            ).all():
                bool_res_not = False
        self.assertTrue(bool_res_not)

    def test_faces_vertices_connectivity_matrix(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        p = MeshFileContainer(geof_file_path)
        faces_vertices_connectivity_matrix = [
            np.array([0, 1], dtype=np.uint8),
            np.array([1, 3], dtype=np.uint8),
            np.array([0, 3], dtype=np.uint8),
            np.array([1, 4], dtype=np.uint8),
            np.array([3, 4], dtype=np.uint8),
            np.array([4, 6], dtype=np.uint8),
            np.array([3, 6], dtype=np.uint8),
            np.array([4, 7], dtype=np.uint8),
            np.array([6, 7], dtype=np.uint8),
            np.array([1, 2], dtype=np.uint8),
            np.array([2, 4], dtype=np.uint8),
            np.array([2, 5], dtype=np.uint8),
            np.array([4, 5], dtype=np.uint8),
            np.array([5, 7], dtype=np.uint8),
            np.array([5, 8], dtype=np.uint8),
            np.array([7, 8], dtype=np.uint8),
        ]
        bool_res_not = True
        for i in range(len(faces_vertices_connectivity_matrix)):
            if (
                p.faces_vertices_connectivity_matrix[i]
                - faces_vertices_connectivity_matrix[i]
            ).all():
                bool_res_not = False
        self.assertTrue(bool_res_not)

    def test_cells_faces_connectivity_matrix(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        p = MeshFileContainer(geof_file_path)
        cells_faces_connectivity_matrix = [
            np.array([1.0, 2.0, 3.0]),
            np.array([4.0, 5.0, 1.0]),
            np.array([4.0, 6.0, 7.0]),
            np.array([8.0, 9.0, 5.0]),
            np.array([10.0, 11.0, 3.0]),
            np.array([12.0, 13.0, 10.0]),
            np.array([12.0, 14.0, 7.0]),
            np.array([15.0, 16.0, 13.0]),
        ]
        bool_res_not = True
        for i in range(len(cells_faces_connectivity_matrix)):
            if (
                p.cells_faces_connectivity_matrix[i]
                - cells_faces_connectivity_matrix[i]
            ).all():
                bool_res_not = False
        self.assertTrue(bool_res_not)

    def test_cells_connectivity_matrix(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        p = MeshFileContainer(geof_file_path)
        cells_connectivity_matrix = [
            np.array([[0, 1], [1, 2], [2, 0]]),
            np.array([[0, 1], [1, 2], [2, 0]]),
            np.array([[0, 1], [1, 2], [2, 0]]),
            np.array([[0, 1], [1, 2], [2, 0]]),
            np.array([[0, 1], [1, 2], [2, 0]]),
            np.array([[0, 1], [1, 2], [2, 0]]),
            np.array([[0, 1], [1, 2], [2, 0]]),
            np.array([[0, 1], [1, 2], [2, 0]]),
        ]
        bool_res_not = True
        for i in range(len(cells_connectivity_matrix)):
            if (p.cells_connectivity_matrix[i] - cells_connectivity_matrix[i]).all():
                bool_res_not = False
        self.assertTrue(bool_res_not)

    def test_vertices_boundaries_connectivity(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        p = MeshFileContainer(geof_file_path)
        vertices_boundaries_connectivity = {
            "BOTTOM": np.array([0, 1, 2], dtype=np.int8),
            "TOP": np.array([6, 7, 8], dtype=np.int8),
            "LEFT": np.array([0, 3, 6], dtype=np.int8),
            "RIGHT": np.array([2, 5, 8], dtype=np.int8),
        }
        bool_res_not = True
        for boundary_name in vertices_boundaries_connectivity.keys():
            if (
                np.array(p.vertices_boundaries_connectivity[boundary_name])
                - np.array(vertices_boundaries_connectivity[boundary_name])
            ).all():
                bool_res_not = False
        self.assertTrue(bool_res_not)

    def test_faces_boundaries_connectivity(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        p = MeshFileContainer(geof_file_path)
        faces_boundaries_connectivity = {
            "BOTTOM": [0, 9],
            "TOP": [8, 15],
            "LEFT": [2, 6],
            "RIGHT": [11, 14],
        }
        bool_res_not = True
        for boundary_name in faces_boundaries_connectivity.keys():
            if (
                np.array(p.faces_boundaries_connectivity[boundary_name])
                - np.array(faces_boundaries_connectivity[boundary_name])
            ).all():
                bool_res_not = False
        self.assertTrue(bool_res_not)

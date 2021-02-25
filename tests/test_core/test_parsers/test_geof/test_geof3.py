from unittest import TestCase

from data.test_data.test_parsers.config import test_parsers_folder_path
from pythhon.parameters import *
import pythhon.mesh.parsers.geof as gf2


class Test(TestCase):
    def test_vertices(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        vertices_check = np.array(
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
            ],
            dtype=real,
        ).T
        vertices = gf2.get_vertices(geof_file_path)
        print(vertices)
        res = not (vertices - vertices_check).all()
        self.assertTrue(res)

    def test_vertices_boundaries_connectivity(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        vertices_boundaries_connectivity_check = {
            "BOTTOM": [0, 1, 2],
            "TOP": [6, 7, 8],
            "LEFT": [0, 3, 6],
            "RIGHT": [2, 5, 8],
        }
        vertices_boundaries_connectivity = gf2.get_vertices_boundaries_connectivity(geof_file_path)
        self.assertTrue(vertices_boundaries_connectivity == vertices_boundaries_connectivity_check)

    def test_get_cells_data(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        a = [
            [0, 1, 3],
            [1, 4, 3],
            [3, 4, 6],
            [4, 7, 6],
            [1, 2, 4],
            [2, 5, 4],
            [4, 5, 7],
            [5, 8, 7],
        ]
        b, _, _ = gf2.get_cells_data(geof_file_path)
        self.assertTrue(b == a)

    def test_get_cells_labels(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        # a = [
        #     "c2d3",
        #     "c2d3",
        #     "c2d3",
        #     "c2d3",
        #     "c2d3",
        #     "c2d3",
        #     "c2d3",
        #     "c2d3",
        # ]
        a = [
            ShapeType.TRIANGLE,
            ShapeType.TRIANGLE,
            ShapeType.TRIANGLE,
            ShapeType.TRIANGLE,
            ShapeType.TRIANGLE,
            ShapeType.TRIANGLE,
            ShapeType.TRIANGLE,
            ShapeType.TRIANGLE,
        ]
        _, b, _ = gf2.get_cells_data(geof_file_path)
        self.assertTrue(b == a)

    def test_vertices_faces_collections(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        a, _, h = gf2.get_cells_data(geof_file_path)
        c, d, _ = gf2.get_faces_data(a, h)
        vertices_faces_connection = [
            [0, 1],
            [1, 3],
            [0, 3],
            [1, 4],
            [3, 4],
            [4, 6],
            [3, 6],
            [4, 7],
            [6, 7],
            [1, 2],
            [2, 4],
            [2, 5],
            [4, 5],
            [5, 7],
            [5, 8],
            [7, 8],
        ]
        print(d)
        self.assertTrue(c == vertices_faces_connection)

    def test_cell_faces_collections(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        a, _, b = gf2.get_cells_data(geof_file_path)
        _, d, _ = gf2.get_faces_data(a, b)
        cells_faces_connection = [
            [0, 1, 2],
            [3, 4, 1],
            [4, 5, 6],
            [7, 8, 5],
            [9, 10, 3],
            [11, 12, 10],
            [12, 13, 7],
            [14, 15, 13],
        ]
        self.assertTrue(d == cells_faces_connection)

    def test_faces_labels(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        a, _, b = gf2.get_cells_data(geof_file_path)
        _, _, d = gf2.get_faces_data(a, b)
        # faces_labels = [
        #     "c1d2",
        #     "c1d2",
        #     "c1d2",
        #     "c1d2",
        #     "c1d2",
        #     "c1d2",
        #     "c1d2",
        #     "c1d2",
        #     "c1d2",
        #     "c1d2",
        #     "c1d2",
        #     "c1d2",
        #     "c1d2",
        #     "c1d2",
        #     "c1d2",
        #     "c1d2",
        # ]
        faces_labels = [
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
            ShapeType.SEGMENT,
        ]
        self.assertTrue(d == faces_labels)

    def test_get_faces_boundaries_connectivity(self):
        geof_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        a, _, b = gf2.get_cells_data(geof_file_path)
        c, d, _ = gf2.get_faces_data(a, b)
        v = gf2.get_vertices_boundaries_connectivity(geof_file_path)
        e = gf2.get_faces_boundaries_connectivity(v, c)
        faces_boundaries_connectivity = {
            "BOTTOM": [0, 9],
            "TOP": [8, 15],
            "LEFT": [2, 6],
            "RIGHT": [11, 14],
        }
        self.assertTrue(e == faces_boundaries_connectivity)

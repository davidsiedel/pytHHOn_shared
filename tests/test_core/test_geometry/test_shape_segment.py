from unittest import TestCase

from parameterized import parameterized
from scipy import integrate

from pythhon.parameters import *
from pythhon.geometry.shapes.segment import *
from pythhon.geometry.shape import Shape


class Test(TestCase):
    def test_segment_0(self):
        v0 = np.array([0.0, 0.0])
        v1 = np.array([0.0, 1.1])
        segment_vertices = np.zeros((2, 2))
        segment_vertices[:, 0] = v0
        segment_vertices[:, 1] = v1
        segment = Shape(ShapeType.SEGMENT, segment_vertices)
        segment_centroid_check = np.array([(v0[0] + v1[0]) / 2.0, (v0[1] + v1[1]) / 2.0])
        segment_volume_check = np.sqrt((v0[0] - v1[0]) ** 2 + (v0[1] - v1[1]) ** 2)
        segment_diameter_check = np.sqrt((v0[0] - v1[0]) ** 2 + (v0[1] - v1[1]) ** 2)
        np.testing.assert_allclose(segment.centroid, segment_centroid_check, rtol=1e-010)
        np.testing.assert_allclose(segment.volume, segment_volume_check, rtol=1e-010)
        np.testing.assert_allclose(segment.diameter, segment_diameter_check, rtol=1e-010)
        self.assertTrue(True)

    def test_segment_1(self):
        v0 = np.array([1.2, 3.4])
        v1 = np.array([-3.6, -3.1])
        segment_vertices = np.zeros((2, 2))
        segment_vertices[:, 0] = v0
        segment_vertices[:, 1] = v1
        segment = Shape(ShapeType.SEGMENT, segment_vertices)
        segment_centroid_check = np.array([(v0[0] + v1[0]) / 2.0, (v0[1] + v1[1]) / 2.0])
        segment_volume_check = np.sqrt((v0[0] - v1[0]) ** 2 + (v0[1] - v1[1]) ** 2)
        segment_diameter_check = np.sqrt((v0[0] - v1[0]) ** 2 + (v0[1] - v1[1]) ** 2)
        np.testing.assert_allclose(segment.centroid, segment_centroid_check, rtol=1e-010)
        np.testing.assert_allclose(segment.volume, segment_volume_check, rtol=1e-010)
        np.testing.assert_allclose(segment.diameter, segment_diameter_check, rtol=1e-010)
        self.assertTrue(True)

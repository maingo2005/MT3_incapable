
#include <Novice.h>
#include <imgui.h>
#include <cmath>
#include <numbers>
#include <algorithm>

struct Vector3 final {
	float x;
	float y;
	float z;
};

struct Matrix4x4 final {
	float m[4][4];
};

struct sphere_s {
	Vector3 center;
	float radius;
};

struct line_s {
	Vector3 origin;
	Vector3 diff;
};

struct ray_s {
	Vector3 origin;
	Vector3 diff;
};

struct segment_s {
	Vector3 origin;
	Vector3 diff;
};

struct plane_s {
	Vector3 normal;
	float distance;
};

struct triangle_s {
	Vector3 vertices[3];
};

struct aabb_s {
	Vector3 min;
	Vector3 max;
};

Vector3 Add(const Vector3& v1, const Vector3& v2) {
	Vector3 result = { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
	return result;
}

Vector3 Subtract(const Vector3& v1, const Vector3& v2) {
	Vector3 result = { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
	return result;
}
Vector3 Multiply(float scalar, const Vector3& v) {
	Vector3 result = { scalar * v.x, scalar * v.y, scalar * v.z };
	return result;
}

float Dot(const Vector3& v1, const Vector3& v2) {
	float result = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	return result;
}
float Length(const Vector3& v) {
	float result = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	return result;
}

Vector3 Normalize(const Vector3& v) {
	float length = std::sqrt(Dot(v, v));
	return Multiply(1.0f / length, v);
}

Vector3 Scale(const Vector3& v, float scale) {
	return { v.x * scale, v.y * scale, v.z * scale };
}

Matrix4x4 Add(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result = {};
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			result.m[i][j] = m1.m[i][j] + m2.m[i][j];
		}
	}
	return result;
}

Matrix4x4 Subtract(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result = {};
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			result.m[i][j] = m1.m[i][j] - m2.m[i][j];
		}
	}
	return result;
}

Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result = {};
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			result.m[i][j] = 0.0f;
			for (int k = 0; k < 4; k++) {
				result.m[i][j] += m1.m[i][k] * m2.m[k][j];
			}
		}
	}
	return result;
}

Matrix4x4 Inverse(const Matrix4x4& m) {

	Matrix4x4 result = {};

	float det =
		m.m[0][0] * (m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[1][2] * m.m[2][3] * m.m[3][1] +
			m.m[1][3] * m.m[2][1] * m.m[3][2] - m.m[1][3] * m.m[2][2] * m.m[3][1] -
			m.m[1][2] * m.m[2][1] * m.m[3][3] - m.m[1][1] * m.m[2][3] * m.m[3][2]) -
		m.m[0][1] * (m.m[1][0] * m.m[2][2] * m.m[3][3] + m.m[1][2] * m.m[2][3] * m.m[3][0] +
			m.m[1][3] * m.m[2][0] * m.m[3][2] - m.m[1][3] * m.m[2][2] * m.m[3][0] -
			m.m[1][2] * m.m[2][0] * m.m[3][3] - m.m[1][0] * m.m[2][3] * m.m[3][2]) +
		m.m[0][2] * (m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[1][1] * m.m[2][3] * m.m[3][0] +
			m.m[1][3] * m.m[2][0] * m.m[3][1] - m.m[1][3] * m.m[2][1] * m.m[3][0] -
			m.m[1][1] * m.m[2][0] * m.m[3][3] - m.m[1][0] * m.m[2][3] * m.m[3][1]) -
		m.m[0][3] * (m.m[1][0] * m.m[2][1] * m.m[3][2] + m.m[1][1] * m.m[2][2] * m.m[3][0] +
			m.m[1][2] * m.m[2][0] * m.m[3][1] - m.m[1][2] * m.m[2][1] * m.m[3][0] -
			m.m[1][1] * m.m[2][0] * m.m[3][2] - m.m[1][0] * m.m[2][2] * m.m[3][1]);

	float recpDeterminant = 1.0f / det;

	result.m[0][0] = (m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[1][2] * m.m[2][3] * m.m[3][1] +
		m.m[1][3] * m.m[2][1] * m.m[3][2] - m.m[1][3] * m.m[2][2] * m.m[3][1] -
		m.m[1][2] * m.m[2][1] * m.m[3][3] - m.m[1][1] * m.m[2][3] * m.m[3][2]) *
		recpDeterminant;
	result.m[0][1] = (-m.m[0][1] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[2][3] * m.m[3][1] -
		m.m[0][3] * m.m[2][1] * m.m[3][2] + m.m[0][3] * m.m[2][2] * m.m[3][1] +
		m.m[0][2] * m.m[2][1] * m.m[3][3] + m.m[0][1] * m.m[2][3] * m.m[3][2]) *
		recpDeterminant;
	result.m[0][2] = (m.m[0][1] * m.m[1][2] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[3][1] +
		m.m[0][3] * m.m[1][1] * m.m[3][2] - m.m[0][3] * m.m[1][2] * m.m[3][1] -
		m.m[0][2] * m.m[1][1] * m.m[3][3] - m.m[0][1] * m.m[1][3] * m.m[3][2]) *
		recpDeterminant;
	result.m[0][3] = (-m.m[0][1] * m.m[1][2] * m.m[2][3] - m.m[0][2] * m.m[1][3] * m.m[2][1] -
		m.m[0][3] * m.m[1][1] * m.m[2][2] + m.m[0][3] * m.m[1][2] * m.m[2][1] +
		m.m[0][2] * m.m[1][1] * m.m[2][3] + m.m[0][1] * m.m[1][3] * m.m[2][2]) *
		recpDeterminant;

	result.m[1][0] = (-m.m[1][0] * m.m[2][2] * m.m[3][3] - m.m[1][2] * m.m[2][3] * m.m[3][0] -
		m.m[1][3] * m.m[2][0] * m.m[3][2] + m.m[1][3] * m.m[2][2] * m.m[3][0] +
		m.m[1][2] * m.m[2][0] * m.m[3][3] + m.m[1][0] * m.m[2][3] * m.m[3][2]) *
		recpDeterminant;
	result.m[1][1] = (m.m[0][0] * m.m[2][2] * m.m[3][3] + m.m[0][2] * m.m[2][3] * m.m[3][0] +
		m.m[0][3] * m.m[2][0] * m.m[3][2] - m.m[0][3] * m.m[2][2] * m.m[3][0] -
		m.m[0][2] * m.m[2][0] * m.m[3][3] - m.m[0][0] * m.m[2][3] * m.m[3][2]) *
		recpDeterminant;
	result.m[1][2] = (-m.m[0][0] * m.m[1][2] * m.m[3][3] - m.m[0][2] * m.m[1][3] * m.m[3][0] -
		m.m[0][3] * m.m[1][0] * m.m[3][2] + m.m[0][3] * m.m[1][2] * m.m[3][0] +
		m.m[0][2] * m.m[1][0] * m.m[3][3] + m.m[0][0] * m.m[1][3] * m.m[3][2]) *
		recpDeterminant;
	result.m[1][3] = (m.m[0][0] * m.m[1][2] * m.m[2][3] + m.m[0][2] * m.m[1][3] * m.m[2][0] +
		m.m[0][3] * m.m[1][0] * m.m[2][2] - m.m[0][3] * m.m[1][2] * m.m[2][0] -
		m.m[0][2] * m.m[1][0] * m.m[2][3] - m.m[0][0] * m.m[1][3] * m.m[2][2]) *
		recpDeterminant;

	result.m[2][0] = (m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[1][1] * m.m[2][3] * m.m[3][0] +
		m.m[1][3] * m.m[2][0] * m.m[3][1] - m.m[1][3] * m.m[2][1] * m.m[3][0] -
		m.m[1][1] * m.m[2][0] * m.m[3][3] - m.m[1][0] * m.m[2][3] * m.m[3][1]) *
		recpDeterminant;
	result.m[2][1] = (-m.m[0][0] * m.m[2][1] * m.m[3][3] - m.m[0][1] * m.m[2][3] * m.m[3][0] -
		m.m[0][3] * m.m[2][0] * m.m[3][1] + m.m[0][3] * m.m[2][1] * m.m[3][0] +
		m.m[0][1] * m.m[2][0] * m.m[3][3] + m.m[0][0] * m.m[2][3] * m.m[3][1]) *
		recpDeterminant;
	result.m[2][2] = (m.m[0][0] * m.m[1][1] * m.m[3][3] + m.m[0][1] * m.m[1][3] * m.m[3][0] +
		m.m[0][3] * m.m[1][0] * m.m[3][1] - m.m[0][3] * m.m[1][1] * m.m[3][0] -
		m.m[0][1] * m.m[1][0] * m.m[3][3] - m.m[0][0] * m.m[1][3] * m.m[3][1]) *
		recpDeterminant;
	result.m[2][3] = (-m.m[0][0] * m.m[1][1] * m.m[2][3] - m.m[0][1] * m.m[1][3] * m.m[2][0] -
		m.m[0][3] * m.m[1][0] * m.m[2][1] + m.m[0][3] * m.m[1][1] * m.m[2][0] +
		m.m[0][1] * m.m[1][0] * m.m[2][3] + m.m[0][0] * m.m[1][3] * m.m[2][1]) *
		recpDeterminant;

	result.m[3][0] = (-m.m[1][0] * m.m[2][1] * m.m[3][2] - m.m[1][1] * m.m[2][2] * m.m[3][0] -
		m.m[1][2] * m.m[2][0] * m.m[3][1] + m.m[1][2] * m.m[2][1] * m.m[3][0] +
		m.m[1][1] * m.m[2][0] * m.m[3][2] + m.m[1][0] * m.m[2][2] * m.m[3][1]) *
		recpDeterminant;
	result.m[3][1] = (m.m[0][0] * m.m[2][1] * m.m[3][2] + m.m[0][1] * m.m[2][2] * m.m[3][0] +
		m.m[0][2] * m.m[2][0] * m.m[3][1] - m.m[0][2] * m.m[2][1] * m.m[3][0] -
		m.m[0][1] * m.m[2][0] * m.m[3][2] - m.m[0][0] * m.m[2][2] * m.m[3][1]) *
		recpDeterminant;
	result.m[3][2] = (-m.m[0][0] * m.m[1][1] * m.m[3][2] - m.m[0][1] * m.m[1][2] * m.m[3][0] -
		m.m[0][2] * m.m[1][0] * m.m[3][1] + m.m[0][2] * m.m[1][1] * m.m[3][0] +
		m.m[0][1] * m.m[1][0] * m.m[3][2] + m.m[0][0] * m.m[1][2] * m.m[3][1]) *
		recpDeterminant;
	result.m[3][3] = (m.m[0][0] * m.m[1][1] * m.m[2][2] + m.m[0][1] * m.m[1][2] * m.m[2][0] +
		m.m[0][2] * m.m[1][0] * m.m[2][1] - m.m[0][2] * m.m[1][1] * m.m[2][0] -
		m.m[0][1] * m.m[1][0] * m.m[2][2] - m.m[0][0] * m.m[1][2] * m.m[2][1]) *
		recpDeterminant;

	return result;
}

Matrix4x4 Transpose(const Matrix4x4& m) {
	Matrix4x4 result = {};
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			result.m[i][j] = m.m[j][i];
		}
	}
	return result;
}

Matrix4x4 MakeIdentity4x4() {
	Matrix4x4 identity = {};
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			identity.m[i][j] = (i == j) ? 1.0f : 0.0f;
		}
	}
	return identity;
}

Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {
	float cotHalfFovV = 1.0f / std::tan(fovY / 2.0f);
	return {
		(cotHalfFovV / aspectRatio),
		0.0f,
		0.0f,
		0.0f,
		0.0f,
		cotHalfFovV,
		0.0f,
		0.0f,
		0.0f,
		0.0f,
		farClip / (farClip - nearClip),
		1.0f,
		0.0f,
		0.0f,
		-(nearClip * farClip) / (farClip - nearClip),
		0.0f };
}

Matrix4x4 MakeRotateXMatrix(float radian) {
	float cosTheta = std::cos(radian);
	float sinTheta = std::sin(radian);
	return { 1.0f, 0.0f,      0.0f,     0.0f,
			0.0f, cosTheta, sinTheta, 0.0f,
			0.0f, -sinTheta, cosTheta, 0.0f,
			0.0f, 0.0f,     0.0f,     1.0f };
}

Matrix4x4 MakeRotateYMatrix(float radian) {
	float cosTheta = std::cos(radian);
	float sinTheta = std::sin(radian);
	return { cosTheta, 0.0f, -sinTheta, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			sinTheta, 0.0f, cosTheta,  0.0f,
			0.0f, 0.0f, 0.0f, 1.0f };
}

Matrix4x4 MakeRotateZMatrix(float radian) {
	float cosTheta = std::cos(radian);
	float sinTheta = std::sin(radian);
	return { cosTheta, sinTheta, 0.0f, 0.0f,
			-sinTheta, cosTheta, 0.0f, 0.0f,
			0.0f,     0.0f,     1.0f, 0.0f,
			0.0f,      0.0f,     0.0f, 1.0f };
}

Matrix4x4 MakeViewportMatrix(
	float left, float top, float width, float height, float minDepth, float maxDepth) {
	return {
		width / 2.0f,
		0.0f,
		0.0f,
		0.0f,
		0.0f,
		-height / 2.0f,
		0.0f,
		0.0f,
		0.0f,
		0.0f,
		maxDepth - minDepth,
		0.0f,
		left + width / 2.0f,
		top + height / 2.0f,
		minDepth,
		1.0f,
	};
}


Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rot, const Vector3& translate) {

	Matrix4x4 ScallMat, RotateMat, RotateMatX, RotateMatY, RotateMatZ, TranslateMat, returnMat;

	// スケール行列作成
	ScallMat = { scale.x, 0, 0, 0, 0, scale.y, 0, 0, 0, 0, scale.z, 0, 0, 0, 0, 1 };

	// XYZ回転行列作成
	RotateMatX = { 1, 0, 0, 0, 0, cosf(rot.x), sinf(rot.x), 0, 0, -sinf(rot.x), cosf(rot.x),
				  0, 0, 0, 0, 1 };

	RotateMatY = { cosf(rot.y), 0, -sinf(rot.y), 0, 0, 1, 0, 0,
				  sinf(rot.y), 0, cosf(rot.y),  0, 0, 0, 0, 1 };

	RotateMatZ = { cosf(rot.z), sinf(rot.z), 0, 0, -sinf(rot.z), cosf(rot.z), 0, 0, 0, 0, 1, 0,
				  0,           0,           0, 1 };

	// XYZ回転行列の合成(Z*X*Y)
	RotateMat = Multiply(RotateMatZ, RotateMatX);
	// ↑の結果＊Y軸回転
	RotateMat = Multiply(RotateMat, RotateMatY);

	// 平行移動行列作成
	TranslateMat = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, translate.x, translate.y, translate.z, 1 };

	// スケール＊回転＊平行移動をワールド変換行列に
	returnMat = Multiply(ScallMat, RotateMat);
	returnMat = Multiply(returnMat, TranslateMat);

	return returnMat;
}

Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {
	Vector3 result;
	// w=1がデカルト座標系であるので(x,y,1)のベクトルとして
	// matrixとの積をとる
	result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] +
		1.0f * matrix.m[3][0];

	result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] +
		1.0f * matrix.m[3][1];

	result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] +
		1.0f * matrix.m[3][2];

	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] +
		1.0f * matrix.m[3][3];

	// w=1がデカルト座標系であるので、w除算することで同次座標をデカルト座標に戻す
	result.x /= w;
	result.y /= w;
	result.z /= w;

	return result;
}

void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	const float kGridHalfWidth = 2.0f;
	const uint32_t kSubdivision = 10;
	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision);
	int color = 0xAAAAAAFF;
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; xIndex++) {
		float x = -kGridHalfWidth + xIndex * kGridEvery;

		Vector3 start(x, 0, -kGridHalfWidth);
		Vector3 end(x, 0, kGridHalfWidth);

		start = Transform(start, viewProjectionMatrix);
		end = Transform(end, viewProjectionMatrix);

		start = Transform(start, viewportMatrix);
		end = Transform(end, viewportMatrix);

		if (xIndex == kSubdivision / 2) {
			color = BLACK;
		}
		else {
			color = 0xAAAAAAFF;
		}
		Novice::DrawLine((int)start.x, (int)start.y, (int)end.x, (int)end.y, color);
	}

	for (uint32_t zIndex = 0; zIndex <= kSubdivision; zIndex++) {
		float z = -kGridHalfWidth + zIndex * kGridEvery;

		Vector3 start(-kGridHalfWidth, 0, z);
		Vector3 end(kGridHalfWidth, 0, z);

		start = Transform(start, viewProjectionMatrix);
		end = Transform(end, viewProjectionMatrix);

		start = Transform(start, viewportMatrix);
		end = Transform(end, viewportMatrix);

		if (zIndex == kSubdivision / 2) {
			color = BLACK;
		}
		else {
			color = 0xAAAAAAFF;
		}
		Novice::DrawLine((int)start.x, (int)start.y, (int)end.x, (int)end.y, color);
	}
}

void DrawSphere(const sphere_s& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	const uint32_t kSubdivision = 20;
	const float kLonEvery = std::numbers::pi_v<float> *2.0f / kSubdivision;
	const float kLatEvery = std::numbers::pi_v<float> / kSubdivision;
	for (uint32_t latIndex = 0; latIndex < kSubdivision; latIndex++) {
		float lat = -std::numbers::pi_v<float> / 2.0f + kLatEvery * latIndex;
		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; lonIndex++) {
			float lon = lonIndex * kLonEvery;
			Vector3 a, b, c;
			a.x = sphere.center.x + sphere.radius * cos(lat) * cos(lon);
			a.y = sphere.center.y + sphere.radius * sin(lat);
			a.z = sphere.center.z + sphere.radius * cos(lat) * sin(lon);
			b.x = sphere.center.x + sphere.radius * cos(lat) * cos(lon + kLonEvery);
			b.y = sphere.center.y + sphere.radius * sin(lat);
			b.z = sphere.center.z + sphere.radius * cos(lat) * sin(lon + kLonEvery);
			c.x = sphere.center.x + sphere.radius * cos(lat + kLatEvery) * cos(lon);
			c.y = sphere.center.y + sphere.radius * sin(lat + kLatEvery);
			c.z = sphere.center.z + sphere.radius * cos(lat + kLatEvery) * sin(lon);

			a = Transform(a, viewProjectionMatrix);
			b = Transform(b, viewProjectionMatrix);
			c = Transform(c, viewProjectionMatrix);

			a = Transform(a, viewportMatrix);
			b = Transform(b, viewportMatrix);
			c = Transform(c, viewportMatrix);
			Novice::DrawLine((int)a.x, (int)a.y, (int)b.x, (int)b.y, color);
			Novice::DrawLine((int)a.x, (int)a.y, (int)c.x, (int)c.y, color);
		}
	}
}

Vector3 Cross(const Vector3& v1, const Vector3& v2) {
	return Vector3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

Vector3 Project(const Vector3& v1, const Vector3& v2) {
	float dotProduct = Dot(v1, v2);
	float lengthSq = Dot(v2, v2);
	return Multiply(dotProduct / lengthSq, v2);
}

Vector3 ClosestPoint(const Vector3& point, const segment_s& segment) {
	Vector3 segmentToPoint = Subtract(point, segment.origin);
	Vector3 projection = Project(segmentToPoint, segment.diff);

	float projectionLengthSq = Dot(projection, projection);
	float segmentLengthSq = Dot(segment.diff, segment.diff);

	if (projectionLengthSq > segmentLengthSq) {
		return Add(segment.origin, segment.diff);
	}
	else if (projectionLengthSq < 0) {
		return segment.origin;
	}
	else {
		return Add(segment.origin, projection);
	}
}

Vector3 Perpendicular(const Vector3& vector) {
	if (vector.x != 0.0f || vector.y != 0.0f) {
		return { -vector.y, vector.x, 0.0f };
	}
	return { 0.0f, -vector.z, vector.y };
}

void DrawPlane(const plane_s& plane, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	Vector3 center = Multiply(plane.distance, plane.normal);
	Vector3 perpendiculars[4];
	perpendiculars[0] = Normalize(Perpendicular(plane.normal));
	perpendiculars[1] = Cross(plane.normal, perpendiculars[0]);
	perpendiculars[2] = { -perpendiculars[0].x, -perpendiculars[0].y, -perpendiculars[0].z };
	perpendiculars[3] = { -perpendiculars[1].x, -perpendiculars[1].y, -perpendiculars[1].z };
	Vector3 points[4];
	for (int32_t index = 0; index < 4; index++) {
		Vector3 extend = Multiply(2.0f, perpendiculars[index]);
		Vector3 point = Add(center, extend);
		points[index] = Transform(Transform(point, viewProjectionMatrix), viewportMatrix);
	}
	for (int32_t i = 0; i < 4; i++) {
		int32_t next = (i + 1) % 4;
		Novice::DrawLine((int)points[i].x, (int)points[i].y, (int)points[next].x, (int)points[next].y, color);
	}
}

void DrawSegment(const segment_s& segment, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t segmentColor) {
	Vector3 start = Transform(Transform(segment.origin, viewProjectionMatrix), viewportMatrix);
	Vector3 end = Transform(Transform(Add(segment.origin, segment.diff), viewProjectionMatrix), viewportMatrix);
	Novice::DrawLine((int)start.x, (int)start.y, (int)end.x, (int)end.y, segmentColor);
}

void DrawTriangle(const triangle_s& triangle, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	Vector3 transformedVertices[3];
	for (int i = 0; i < 3; ++i) {
		transformedVertices[i] = Transform(triangle.vertices[i], viewProjectionMatrix);
		transformedVertices[i] = Transform(transformedVertices[i], viewportMatrix);
	}
	for (int i = 0; i < 3; ++i) {
		int next = (i + 1) % 3;
		Novice::DrawLine((int)transformedVertices[i].x, (int)transformedVertices[i].y, (int)transformedVertices[next].x, (int)transformedVertices[next].y, color);
	}
}

void DrawAABB(const aabb_s& aabb, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	Vector3 corners[8] = {
		{aabb.min.x, aabb.min.y, aabb.min.z},
		{aabb.max.x, aabb.min.y, aabb.min.z},
		{aabb.max.x, aabb.max.y, aabb.min.z},
		{aabb.min.x, aabb.max.y, aabb.min.z},
		{aabb.min.x, aabb.min.y, aabb.max.z},
		{aabb.max.x, aabb.min.y, aabb.max.z},
		{aabb.max.x, aabb.max.y, aabb.max.z},
		{aabb.min.x, aabb.max.y, aabb.max.z}
	};

	for (int i = 0; i < 8; ++i) {
		corners[i] = Transform(Transform(corners[i], viewProjectionMatrix), viewportMatrix);
	}

	int index[24] = { 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4, 0, 4, 1, 5, 2, 6, 3, 7 };

	for (int i = 0; i < 24; i += 2) {
		Novice::DrawLine((int)corners[index[i]].x, (int)corners[index[i]].y, (int)corners[index[i + 1]].x, (int)corners[index[i + 1]].y, color);
	}
}

Vector3 operator+(const Vector3& v) { return v; }

Vector3 operator-(const Vector3& v) { return Vector3(-v.x, -v.y, -v.z); }

// 2項演算子オーバーロード
const Vector3 operator+(const Vector3& v1, const Vector3& v2) {
	return Add(v1, v2);
}

const Vector3 operator-(const Vector3& v1, const Vector3& v2) {
	return Subtract(v1, v2);
}

const Vector3 operator*(float s, const Vector3& v) {
	return Multiply(s, v);
}

const Vector3 operator*(const Vector3& v, float s) {
	return s * v;
}

const Vector3 operator/(const Vector3& v, float s) {
	return Multiply(1.0f / s, v);
}

Matrix4x4 operator+(const Matrix4x4& m1, const Matrix4x4& m2) {
	return Add(m1, m2);
}

Matrix4x4 operator-(const Matrix4x4& m1, const Matrix4x4& m2) {
	return Subtract(m1, m2);
}

Matrix4x4 operator*(const Matrix4x4& m1, const Matrix4x4& m2) {
	return Multiply(m1, m2);
}



Vector3& operator+=(Vector3& lhv, const Vector3& rhv) {
	lhv.x += rhv.x;
	lhv.y += rhv.y;
	lhv.z += rhv.z;
	return lhv;
}
Vector3& operator-=(Vector3& lhv, const Vector3& rhv) {
	lhv.x -= rhv.x;
	lhv.y -= rhv.y;
	lhv.z -= rhv.z;
	return lhv;
}
Vector3& operator*=(Vector3& v, float s) {
	v.x *= s;
	v.y *= s;
	v.z *= s;
	return v;
}
Vector3& operator/=(Vector3& v, float s) {
	v.x /= s;
	v.y /= s;
	v.z /= s;
	return v;
}

bool IsCollision(const sphere_s& s1, const sphere_s& s2) {
	float distance = Length(s2.center - s1.center);
	return distance <= s1.radius + s2.radius;

}

bool IsCollision(const sphere_s& sphere, const plane_s& plane) {
	float distance = Dot(plane.normal, sphere.center) - plane.distance;
	return std::abs(distance) <= sphere.radius;
}

bool IsCollision(const segment_s& segment, const plane_s& plane) {
	float dot = Dot(segment.diff, plane.normal);
	if (dot == 0.0f) {
		return false;
	}
	float t = (plane.distance - Dot(segment.origin, plane.normal)) / dot;

	return t >= 0.0f && t <= 1.0f;
}

bool IsCollision(const triangle_s& triangle, const segment_s& segment) {
	Vector3 v0 = triangle.vertices[0];
	Vector3 v1 = triangle.vertices[1];
	Vector3 v2 = triangle.vertices[2];

	Vector3 v01 = Subtract(v1, v0);
	Vector3 v12 = Subtract(v2, v1);
	Vector3 v20 = Subtract(v0, v2);

	Vector3 normal = Cross(v01, v12);

	float dot = Dot(normal, segment.diff);
	if (dot == 0.0f) {
		return false;
	}

	float t = (Dot(normal, v0) - Dot(normal, segment.origin)) / dot;
	if (t < 0.0f || t > 1.0f) {
		return false;
	}

	Vector3 intersectionPoint = Add(segment.origin, Multiply(t, segment.diff));

	Vector3 v1p = Subtract(intersectionPoint, v1);
	Vector3 v2p = Subtract(intersectionPoint, v2);
	Vector3 v0p = Subtract(intersectionPoint, v0);
	Vector3 cross01 = Cross(v01, v1p);
	Vector3 cross12 = Cross(v12, v2p);
	Vector3 cross20 = Cross(v20, v0p);

	return Dot(cross01, normal) >= 0.0f && Dot(cross12, normal) >= 0.0f && Dot(cross20, normal) >= 0.0f;
}

bool IsCollision(const aabb_s& aabb1, const aabb_s& aabb2) {
	return (aabb1.min.x <= aabb2.max.x && aabb1.max.x >= aabb2.min.x) &&
		(aabb1.min.y <= aabb2.max.y && aabb1.max.y >= aabb2.min.y) &&
		(aabb1.min.z <= aabb2.max.z && aabb1.max.z >= aabb2.min.z);
}

bool IsCollision(const aabb_s& aabb, const sphere_s& sphere) {
	Vector3 closestPoint{
		std::clamp(sphere.center.x, aabb.min.x, aabb.max.x),
		std::clamp(sphere.center.y, aabb.min.y, aabb.max.y),
		std::clamp(sphere.center.z, aabb.min.z, aabb.max.z)
	};
	float distance = Length(closestPoint - sphere.center);
	if (distance <= sphere.radius) {
		return true;
	}
	else {
		return false;
	}
}

bool IsCollision(const aabb_s& aabb, const segment_s& segment) {
	float tNearX = (aabb.min.x - segment.origin.x) / segment.diff.x;
	float tFarX = (aabb.max.x - segment.origin.x) / segment.diff.x;

	if (tNearX > tFarX) {
		std::swap(tNearX, tFarX);
	}

	float tNearY = (aabb.min.y - segment.origin.y) / segment.diff.y;
	float tFarY = (aabb.max.y - segment.origin.y) / segment.diff.y;

	if (tNearY > tFarY) {
		std::swap(tNearY, tFarY);
	}

	float tNearZ = (aabb.min.z - segment.origin.z) / segment.diff.z;
	float tFarZ = (aabb.max.z - segment.origin.z) / segment.diff.z;

	if (tNearZ > tFarZ) {
		std::swap(tNearZ, tFarZ);
	}

	float tmin = (std::max)((std::max)(tNearX, tNearY), tNearZ);
	float tmax = (std::min)((std::min)(tFarX, tFarY), tFarZ);

	if (tmin <= tmax && tmax >= 0.0f && tmin <= 1.0f) {
		return true;
	}

	return false;
}

Vector3 Lerp(const Vector3& v1, const Vector3& v2, float t) {
	return { v1.x + t * (v2.x - v1.x),
			v1.y + t * (v2.y - v1.y),
			v1.z + t * (v2.z - v1.z)
	};
}

void DrawBezier(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	const int numSegments = 100;
	Vector3 previousPoint = p0;

	for (int i = 1; i <= numSegments; ++i) {
		float t = static_cast<float>(i) / numSegments;

		Vector3 p0p1 = Lerp(p0, p1, t);
		Vector3 p1p2 = Lerp(p1, p2, t);
		Vector3 currentPoint = Lerp(p0p1, p1p2, t);

		Vector3 start = Transform(Transform(previousPoint, viewProjectionMatrix), viewportMatrix);
		Vector3 end = Transform(Transform(currentPoint, viewProjectionMatrix), viewportMatrix);

		Novice::DrawLine((int)start.x, (int)start.y, (int)end.x, (int)end.y, color);

		previousPoint = currentPoint;
	}
}

const char kWindowTitle[] = "MT3_03_02_Basic";

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	int kWindowWidth = 1280;
	int kWindowHeight = 720;
	Novice::Initialize(kWindowTitle, kWindowWidth, kWindowHeight);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	Vector3 cameraTranslate{ 0.0f, 1.9f, -6.49f };
	Vector3 cameraRotate{ 0.26f, 0.0f, 0.0f };
	//Vector3 rotate{};
	Vector3 translate{};

	sphere_s sphere[2]{
		{.center{0.0f, 0.0f, 0.0f}, .radius{0.5f}},
		{.center{0.9f, 0.0f, 0.0f}, .radius{0.3f}}
	};

	segment_s segment{
		{-0.7f, -0.5f, -0.5f},
		{0.5f,  0.5f,  0.5f}
	};

	plane_s plane{
		.normal = {0.0f, 1.0f, 0.0f},
		.distance = 1.0f
	};

	triangle_s triangle{
		{{-1.0f, 0.0f, 0.0f},
		{0.0f, 1.0f, 0.0f},
		{1.0f, 0.0f, 0.0f}}
	};

	aabb_s aabb1{
		.min{-0.5f, -0.5f, -0.5f},
		.max{0.5f,  0.5f,  0.5f },
	};

	aabb_s aabb2{
		.min{0.2f, 0.2f, 0.2f},
		.max{1.0f, 1.0f, 1.0f},
	};

	Vector3 controlPoint[3] = {
		{-0.8f, 0.58f, 1.0f },
		{1.76f, 1.0f, -0.3f },
		{0.94f, -0.7f, 2.3f }
	};

	Vector3 translates[3] = {
		{0.2f, 1.0f, 0.0f},
		{0.4f, 0.0f, 0.0f},
		{0.3f, 0.0f, 0.0f},
	};

	Vector3 rotates[3] = {
		{0.0f, 0.0f, -6.8f},
		{0.0f, 0.0f, -1.4f},
		{0.0f, 0.0f, 0.0f },
	};

	Vector3 scales[3] = {
		{1.0f, 1.0f, 1.0f},
		{1.0f, 1.0f, 1.0f},
		{1.0f, 1.0f, 1.0f}
	};

	Vector3 point{ -1.5f, 0.6f, 0.6f };
	Vector3 project = Project(Subtract(point, segment.origin), segment.diff);
	Vector3 closestPoint = ClosestPoint(point, segment);
	sphere_s pointSphere{ point, 0.01f };
	sphere_s closestPointSphere{ closestPoint, 0.01f };
	uint32_t sphereColor[2] = { WHITE, WHITE };
	//uint32_t segmentColor = WHITE;
	//uint32_t bezierColor = BLUE;

	Vector3 a{ 0.2f, 1.0f, 0.0f };
	Vector3 b{ 2.4f, 3.1f, 1.2f };

	Vector3 c = a + b;
	Vector3 d = a - b;
	Vector3 e = a * 2.4f;

	Vector3 rotate{ 0.4f, 1.43f, -0.8f };
	Matrix4x4 rotateXMatrix = MakeRotateXMatrix(rotate.x);
	Matrix4x4 rotateYMatrix = MakeRotateYMatrix(rotate.y);
	Matrix4x4 rotateZMatrix = MakeRotateZMatrix(rotate.z);
	Matrix4x4 rotateMatrix = rotateXMatrix * rotateYMatrix * rotateZMatrix;
	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		ImGui::Begin("Window");
		//ImGui::DragFloat3("CameraTranslate", &cameraTranslate.x, 0.01f);
		//ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);
		//ImGui::DragFloat3("Sphere[0]Center", &sphere[0].center.x, 0.01f);
		//ImGui::DragFloat("Sphere[0]Radius", &sphere[0].radius, 0.01f);
		//ImGui::DragFloat3("Sphere[1]Center", &sphere[1].center.x, 0.01f);
		//ImGui::DragFloat("Sphere[1]Radius", &sphere[1].radius, 0.01f);
		//ImGui::DragFloat3("Plane.Normal", &plane.normal.x, 0.01f);
		//ImGui::DragFloat("Plane.Distance", &plane.distance, 0.01f);
		//ImGui::DragFloat3("segment.origin", &segment.origin.x, 0.01f);
		//ImGui::DragFloat3("segment.diff", &segment.diff.x, 0.01f);
		//ImGui::DragFloat3("Traiangle.v0", &triangle.vertices[0].x, 0.01f);
		//ImGui::DragFloat3("Traiangle.v1", &triangle.vertices[1].x, 0.01f);
		//ImGui::DragFloat3("Traiangle.v2", &triangle.vertices[2].x, 0.01f);
		//ImGui::DragFloat3("aabb1.min", &aabb1.min.x, 0.01f);
		//ImGui::DragFloat3("aabb1.max", &aabb1.max.x, 0.01f);
		//ImGui::DragFloat3("aabb2.min", &aabb2.min.x, 0.01f);
		//ImGui::DragFloat3("aabb2.max", &aabb2.max.x, 0.01f);
		//ImGui::InputFloat3("Project", &project.x, "%.3f", ImGuiInputTextFlags_ReadOnly);
		//ImGui::DragFloat3("controlPoint[0]", &controlPoint[0].x, 0.01f);
		//ImGui::DragFloat3("controlPoint[1]", &controlPoint[1].x, 0.01f);
		//ImGui::DragFloat3("controlPoint[2]", &controlPoint[2].x, 0.01f);
		//ImGui::DragFloat3("translates[0]", &translates[0].x, 0.01f);
		//ImGui::DragFloat3("rotates[0]", &rotates[0].x, 0.01f);
		//ImGui::DragFloat3("scales[0]", &scales[0].x, 0.01f);
		//ImGui::DragFloat3("translates[1]", &translates[1].x, 0.01f);
		//ImGui::DragFloat3("rotates[1]", &rotates[1].x, 0.01f);
		//ImGui::DragFloat3("scales[1]", &scales[1].x, 0.01f);
		//ImGui::DragFloat3("translates[2]", &translates[2].x, 0.01f);
		//ImGui::DragFloat3("rotates[2]", &rotates[2].x, 0.01f);
		//ImGui::DragFloat3("scales[2]", &scales[2].x, 0.01f);
		ImGui::Text("c:%f, %f, %f", c.x, c.y, c.z);
		ImGui::Text("d:%f, %f, %f", d.x, d.y, d.z);
		ImGui::Text("e:%f, %f, %f", e.x, e.y, e.z);
		ImGui::Text(
			"matrix:\n%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n",
			rotateMatrix.m[0][0], rotateMatrix.m[0][1], rotateMatrix.m[0][2], rotateMatrix.m[0][3],
			rotateMatrix.m[1][0], rotateMatrix.m[1][1], rotateMatrix.m[1][2], rotateMatrix.m[1][3],
			rotateMatrix.m[2][0], rotateMatrix.m[2][1], rotateMatrix.m[2][2], rotateMatrix.m[2][3],
			rotateMatrix.m[3][0], rotateMatrix.m[3][1], rotateMatrix.m[3][2], rotateMatrix.m[3][3]);
		ImGui::End();

		//aabb1.min.x = (std::min)(aabb1.min.x, aabb1.max.x);
		//aabb1.min.y = (std::min)(aabb1.min.y, aabb1.max.y);
		//aabb1.min.z = (std::min)(aabb1.min.z, aabb1.max.z);
		//aabb1.max.x = (std::max)(aabb1.min.x, aabb1.max.x);
		//aabb1.max.y = (std::max)(aabb1.min.y, aabb1.max.y);
		//aabb1.max.z = (std::max)(aabb1.min.z, aabb1.max.z);

		//aabb2.min.x = (std::min)(aabb2.min.x, aabb2.max.x);
		//aabb2.min.y = (std::min)(aabb2.min.y, aabb2.max.y);
		//aabb2.min.z = (std::min)(aabb2.min.z, aabb2.max.z);
		//aabb2.max.x = (std::max)(aabb2.min.x, aabb2.max.x);
		//aabb2.max.y = (std::max)(aabb2.min.y, aabb2.max.y);
		//aabb2.max.z = (std::max)(aabb2.min.z, aabb2.max.z);

		plane.normal = Normalize(plane.normal);
		/// ↓更新処理ここから
		// カメラの拡大縮小、回転、平行移動ベクトルからカメラ空間の行列を作る
		Matrix4x4 cameraMatrix =
			MakeAffineMatrix({ 1.0f, 1.0f, 1.0f }, cameraRotate, cameraTranslate);

		// ビュー行列：
		// カメラから見た視錐台空間に変換する行列。カメラ行列の逆行列
		Matrix4x4 viewMatrix = Inverse(cameraMatrix);

		// プロジェクション行列：
		// ビュー空間の視錐台を幅と高さが2、奥行きが1の直方体空間に変換する行列
		// この座標系を正規デバイス座標系という
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(
			0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);

		// ビュー行列とプロジェクション行列を掛けて1つの行列に
		// 物体の頂点にこれを掛けると平面に投影された2次元画像の座標になる
		Matrix4x4 viewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);

		// ビューポート行列：
		// 投影面に投影された2次元映像を、スクリーン座標系で表された
		// ウインドウ上の指定領域であるビューポート内に表示する変換行列
		Matrix4x4 viewportMatrix =
			MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

		//Matrix4x4 localMatrix0 = MakeAffineMatrix(scales[0], rotates[0], translates[0]);
		//Matrix4x4 localMatrix1 = MakeAffineMatrix(scales[1], rotates[1], translates[1]);
		//Matrix4x4 localMatrix2 = MakeAffineMatrix(scales[2], rotates[2], translates[2]);
		//
		//Matrix4x4 worldMatrix0 = localMatrix0;
		//Matrix4x4 worldMatrix1 = Multiply(localMatrix1, worldMatrix0);
		//Matrix4x4 worldMatrix2 = Multiply(localMatrix2, worldMatrix1);
		//
		//sphere_s sphere0{
		//    .center = {worldMatrix0.m[3][0], worldMatrix0.m[3][1], worldMatrix0.m[3][2]},
		//    .radius = 0.05f,
		//};
		//sphere_s sphere1{
		//    .center = {worldMatrix1.m[3][0], worldMatrix1.m[3][1], worldMatrix1.m[3][2]},
		//    .radius = 0.05f,
		//};
		//sphere_s sphere2{
		//    .center = {worldMatrix2.m[3][0], worldMatrix2.m[3][1], worldMatrix2.m[3][2]},
		//    .radius = 0.05f
		//};

		//segmentColor = IsCollision(triangle, segment) ? RED : WHITE;

		//sphereColor[0] = IsCollision(sphere[0], sphere[1]) ? RED : WHITE;

		//sphereColor[0] = IsCollision(sphere[0], plane) ? RED : WHITE;

		//segmentColor = IsCollision(segment, plane) ? RED : WHITE;

		//int32_t aabbColor = IsCollision(aabb1, aabb2) ? RED : WHITE;

		//int32_t aabbColor = IsCollision(aabb1, sphere[0]) ? RED : WHITE;

		//int32_t aabbColor = IsCollision(aabb1, segment) ? RED : WHITE;

		/// ↑更新処理ここまで



		/// ↓描画処理ここから

		// ワールド空間の頂点×ビュープロジェクション行列×ビューポート行列で
		// スクリーン空間の座標になる
		DrawGrid(viewProjectionMatrix, viewportMatrix);
		//DrawSegment(segment, viewProjectionMatrix, viewportMatrix, WHITE);
		//DrawTriangle(triangle, viewProjectionMatrix, viewportMatrix, WHITE);
		//DrawAABB(aabb1, viewProjectionMatrix, viewportMatrix, aabbColor);
		//DrawAABB(aabb2, viewProjectionMatrix, viewportMatrix, WHITE);
		//DrawSphere(sphere[0], viewProjectionMatrix, viewportMatrix, sphereColor[0]);
		//DrawSphere(sphere[1], viewProjectionMatrix, viewportMatrix, sphereColor[1]);
		//DrawSphere(pointSphere, viewProjectionMatrix, viewportMatrix, RED);
		//DrawSphere(closestPointSphere, viewProjectionMatrix, viewportMatrix, BLACK);
		//DrawPlane(plane, viewProjectionMatrix, viewportMatrix, sphereColor[0]);
		//DrawBezier(controlPoint[0], controlPoint[1], controlPoint[2], viewProjectionMatrix, viewportMatrix, bezierColor);
		//DrawSphere(sphere0, viewProjectionMatrix, viewportMatrix, RED);
		//DrawSphere(sphere1, viewProjectionMatrix, viewportMatrix, GREEN);
		//DrawSphere(sphere2, viewProjectionMatrix, viewportMatrix, BLUE);
		//DrawSegment(segment_s{sphere0.center, sphere1.center - sphere0.center}, viewProjectionMatrix,viewportMatrix, WHITE);
		//DrawSegment(segment_s{sphere1.center, sphere2.center - sphere1.center}, viewProjectionMatrix,viewportMatrix, WHITE);

		/// ↑描画処理ここまで


		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}

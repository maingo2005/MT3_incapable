#include <Novice.h>

//#include "../MT3.h"
#include <imgui.h>

const char kWindowTitle[] = "MT3_04_00_Basic";

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	struct Spring {
		Vector3 anchor;           // アンカーポイント。固定された端の位置
		float naturalLength;      // 自然長
		float stiffness;          // 剛性。バネ定数k
		float dampingCoefficient; // 減衰係数
	};

	struct Ball {
		Vector3 position;     // ボールの位置
		Vector3 velocity;     // ボールの速度
		Vector3 acceleration; // ボールの加速度
		float mass;           // ボールの質量
		float radius;         // ボールの半径
		unsigned int color;   // ボールの色
	};

	// ライブラリの初期化
	const int kWindowWidth = 1280;
	const int kWindowHeight = 720;
	const float kAxesSize = 100.0f;
	Novice::Initialize(kWindowTitle, kWindowWidth, kWindowHeight);

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

	Camera camera{};
	camera.spherical.theta = 0.26f;
	camera.spherical.radius = 7.0f;

	Spring spring{};
	spring.anchor = {0.0f, 0.0f, 0.0f};
	spring.naturalLength = 1.0f;
	spring.stiffness = 100.0f;
	spring.dampingCoefficient = 2.0f;

	Ball ball{};
	ball.position = {1.2f, 0.0f, 0.0f};
	ball.mass = 2.0f;
	ball.radius = 0.05f;
	ball.color = BLUE;
	float deltaTime = 1.0f / 60.0f;
	bool start = false;

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		Update(camera);

		ImGui::Begin("Window");
		if (ImGui::Button("Start")) {
			start = true;
		}
		ImGui::End();

		if (start) {
			Vector3 diff = ball.position - spring.anchor;
			float length = Length(diff);
			if (length != 0.0f) {
				Vector3 direction = Normalize(diff);
				Vector3 restPosition = spring.anchor + direction * spring.naturalLength;
				Vector3 displacement = length * (ball.position - restPosition);
				Vector3 restoringForce = -spring.stiffness * displacement;
				Vector3 dampingForce = -spring.dampingCoefficient * ball.velocity;
				Vector3 force = restoringForce + dampingForce;
				ball.acceleration = force / ball.mass;
			}
			ball.velocity += ball.acceleration * deltaTime;
			ball.position += ball.velocity * deltaTime;
		}

		Matrix4x4 viewMatrix = CalcViewMatrix(camera);
		Matrix4x4 projectionMatrix =
		  MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		Matrix4x4 viewProjectionMatrix = viewMatrix * projectionMatrix;

		// ViewportMatrixを作る
		Matrix4x4 viewportMatrix =
		  MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

		Matrix4x4 viewportMatrixForAxes = MakeViewportMatrix(
		  0.0f, float(kWindowHeight) - kAxesSize, kAxesSize, kAxesSize, 0.0f, 1.0f);

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		// 描画
		DrawGrid(viewProjectionMatrix, viewportMatrix);

		DrawSphere(
		  Sphere{.center = ball.position, .radius = ball.radius}, viewProjectionMatrix,
		  viewportMatrix, ball.color);

		Segment toPosition;
		toPosition.origin = spring.anchor;
		toPosition.diff = ball.position - spring.anchor;
		DrawSegment(toPosition, viewProjectionMatrix, viewportMatrix, WHITE);

		DrawAxis(camera, viewportMatrixForAxes);

		///
		/// ↑描画処理ここまで
		///

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

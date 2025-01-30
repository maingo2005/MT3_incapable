#include <Novice.h>

//#include "../MT3.h"
#include <imgui.h>

const char kWindowTitle[] = "MT3_04_02_Basic";

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

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

	Ball ball{};
	ball.position = {1.2f, 0.0f, 0.0f};
	ball.mass = 2.0f;
	ball.radius = 0.05f;
	ball.color = WHITE;
	float deltaTime = 1.0f / 60.0f;
	bool start = false;

	struct Pendulum {
		Vector3 anchor;            // アンカーポイント。固定された端の位置
		float length;              // 紐の長さ
		float angle;               // 現在の角度
		float angularVelocity;     // 角速度ω
		float angularAcceleration; // 角加速度
	};

	Pendulum pendulum;
	pendulum.anchor = {0.0f, 1.0f, 0.0f};
	pendulum.length = 0.8f;
	pendulum.angle = 0.7f;
	pendulum.angularVelocity = 0.0f;
	pendulum.angularAcceleration = 0.0f;

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
			pendulum.angularAcceleration = -(9.8f / pendulum.length) * std::sin(pendulum.angle);
			pendulum.angularVelocity += pendulum.angularAcceleration * deltaTime;
			pendulum.angle += pendulum.angularVelocity * deltaTime;
		}
		ball.position.x = pendulum.anchor.x + std::sin(pendulum.angle) * pendulum.length;
		ball.position.y = pendulum.anchor.y - std::cos(pendulum.angle) * pendulum.length;
		ball.position.z = pendulum.anchor.z;

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

		DrawSegment(
		  Segment{pendulum.anchor, ball.position - pendulum.anchor}, viewProjectionMatrix,
		  viewportMatrix, WHITE);

		DrawSphere(
		  Sphere{.center = ball.position, .radius = ball.radius}, viewProjectionMatrix,
		  viewportMatrix, ball.color);

		// DrawSphere(
		//   Sphere{.center = origin, .radius = 0.3f}, viewProjectionMatrix, viewportMatrix, BLUE);

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
